/*
 * Densified exports -> one VQSR-filtered sites-only VCF per VDS, and the whole (genotypes-inclusive) fragments with
 * VQSR's filters and annotations, headless, for a compose to join into one VCF.
 *
 * densified: [plan from planVqsr() in main.nf, header-per-shard manifest, separate-header manifest,
 *            whole fragments' (header-per-shard) manifest]
 * vqsr: settings from main.nf. fragments_per_job, the filter levels, and resources: {name: VCF path}
 *
 * The models are trained on the whole callset, joined by COMPOSE_SITES_ONLY. The SNP tranches are trained, and both
 * modes applied, on chunks of header-per-shard shards. The whole fragments are chunked the same way, and annotated
 * from the matching shards. Everything is keyed by the VDS name, as a run can't tell how many VDSs it will densify
 * until it's planned them.
 */

include { COMPOSE_SITES_ONLY } from '../../../modules/local/compose'
include {
    TRAIN_VQSR_INDEL_MODEL ;
    TRAIN_VQSR_SNP_MODEL ;
    TRAIN_VQSR_SNP_TRANCHES ;
    GATHER_VQSR_SNP_TRANCHES ;
    APPLY_VQSR ;
    CONCAT_VCFS ;
    ANNOTATE_FRAGMENTS
} from '../../../modules/local/vqsr'

// a resource VCF and its index: .tbi for bgzipped, .idx (Tribble) for plain text
def withIndex(path) {
    [file(path), file(path.endsWith('.gz') ? "${path}.tbi" : "${path}.idx")]
}

def stagedResources(vqsr, names) {
    names.collectMany { name -> withIndex(vqsr.resources[name]) }
}

// the shard names in a header-per-shard export's manifest, in genomic order
def manifestNames(manifest) {
    manifest.readLines().collect { line -> line.trim() }.findAll { line -> line }
}

// a header-per-shard export, in chunks: [[VDS name, chunk index, shard VCFs, their indexes]], in genomic order
def chunked(key, dir, names, size) {
    def batches = names.collate(size)
    (0..<batches.size()).collect { index ->
        def shards = batches[index].collect { name -> "${dir}/${name}" }
        [key, index, shards.collect { shard -> file(shard) }, shards.collect { shard -> file("${shard}.tbi") }]
    }
}

workflow VQSR {
    take:
    densified
    vqsr

    main:
    def snp_refs = stagedResources(vqsr, ['hapmap', 'omni', 'one_thousand_genomes', 'dbsnp'])
    def indel_refs = stagedResources(vqsr, ['mills', 'axiom_poly', 'dbsnp'])

    // [VDS name, plan], to put the plan back on outputs keyed by name
    def plans = densified.map { plan, _sites_only, _separate_header, _whole -> [plan.name, plan] }

    COMPOSE_SITES_ONLY(densified.map { plan, _sites_only, separate_header, _whole -> [plan, separate_header] })
    def composed = COMPOSE_SITES_ONLY.out.map { plan -> [plan, file(plan.composed)] }

    TRAIN_VQSR_INDEL_MODEL(composed, indel_refs, vqsr)
    TRAIN_VQSR_SNP_MODEL(composed, snp_refs, vqsr)

    def chunks = densified.flatMap { plan, sites_only, _separate_header, _whole ->
        chunked(plan.name, plan.sites_only, manifestNames(sites_only), vqsr.fragments_per_job)
    }

    // the whole fragments, chunked like the sites-only shards. Both are exported from the dense MT's partitions, so
    // a whole fragment holds the same records as the sites-only shard of the same name
    def whole_chunks = densified.flatMap { plan, sites_only, _separate_header, whole ->
        def names = manifestNames(whole)
        if (names != manifestNames(sites_only)) {
            error("${plan.name}: the whole fragments and sites-only shards are different partitions, so can't be paired")
        }
        chunked(plan.name, plan.whole_fragments, names, vqsr.fragments_per_job)
    }

    TRAIN_VQSR_SNP_TRANCHES(chunks.combine(TRAIN_VQSR_SNP_MODEL.out, by: 0), snp_refs, vqsr)

    // every shard's tranches, gathered once all chunks are trained
    def shard_tranches = TRAIN_VQSR_SNP_TRANCHES.out.tranches
        .groupTuple()
        .map { key, tranches -> [key, tranches.flatten().toSorted { tranche -> tranche.name }] }
    GATHER_VQSR_SNP_TRANCHES(shard_tranches)

    // each chunk meets its own shards' SNP recalibrations, the gathered SNP tranches, and the indel model
    def indel_model = TRAIN_VQSR_INDEL_MODEL.out.map { plan, recal, recal_idx, tranches -> [plan.name, recal, recal_idx, tranches] }
    def to_apply = chunks
        .join(TRAIN_VQSR_SNP_TRANCHES.out.recal, by: [0, 1])
        .combine(GATHER_VQSR_SNP_TRANCHES.out, by: 0)
        .combine(indel_model, by: 0)
    APPLY_VQSR(to_apply, vqsr)

    // the chunks back in genomic order, for a plain concat
    def vqsr_chunks = APPLY_VQSR.out.chunk
        .map { key, index, vcf, tbi -> [key, [index, vcf, tbi]] }
        .groupTuple()
        .join(plans)
        .map { _key, parts, plan ->
            def ordered = parts.toSorted { part -> part[0] }
            [plan, ordered.collect { part -> part[1] }, ordered.collect { part -> part[2] }]
        }
    CONCAT_VCFS(vqsr_chunks)

    // each chunk of whole fragments meets the VQSR VCFs of the same shards
    ANNOTATE_FRAGMENTS(whole_chunks.join(APPLY_VQSR.out.shards, by: [0, 1]))

    // every chunk writes the same header and EOF block: chunk 0's are used
    def headers = ANNOTATE_FRAGMENTS.out.header
        .filter { _key, index, _header, _eof -> index == 0 }
        .map { key, _index, header, eof -> [key, header, eof] }

    // the fragments in genomic order: chunks by index, and within one, by partition name (part-NNNNN, zero-padded).
    // A chunk of one fragment emits a lone path, hence the flatten
    def annotated_fragments = ANNOTATE_FRAGMENTS.out.fragments
        .map { key, index, fragments -> [key, [index, [fragments].flatten().toSorted { fragment -> fragment.name }]] }
        .groupTuple()
        .join(headers)
        .join(plans)
        .map { _key, parts, header, eof, plan ->
            [plan, header, parts.toSorted { part -> part[0] }.collectMany { part -> part[1] }, eof]
        }

    emit:
    // [plan, VQSR sites-only VCF, its index]
    vcf = CONCAT_VCFS.out
    // [plan, header, headless whole fragments in genomic order, BGZF EOF block]: composed in that order, one VCF
    annotated = annotated_fragments
}
