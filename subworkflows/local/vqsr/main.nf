/*
 * Densified exports -> the whole (genotypes-inclusive) fragments with VQSR's filters and annotations, each a VCF with
 * its header, for DATASET_VCFS to split by dataset.
 *
 * densified: [plan from planVqsr() in main.nf, header-per-shard manifest, separate-header manifest,
 *            whole fragments' (header-per-shard) manifest]
 * vqsr: settings from main.nf. fragments_per_job, annotate_fragments_per_job, the filter levels, and resources:
 *       {name: VCF path}
 *
 * The models are trained on the whole callset, joined by COMPOSE_SITES_ONLY. The SNP tranches are trained, and both
 * modes applied, on chunks of fragments_per_job header-per-shard shards. The whole fragments are chunked separately,
 * by annotate_fragments_per_job, and annotated from the matching shards' VQSR VCFs, regrouped to match. Everything
 * is keyed by the VDS name, as a run can't tell how many VDSs it will densify until it's planned them.
 */

include { COMPOSE_SITES_ONLY } from '../../../modules/local/compose'
include {
    TRAIN_VQSR_INDEL_MODEL ;
    TRAIN_VQSR_SNP_MODEL ;
    TRAIN_VQSR_SNP_TRANCHES ;
    GATHER_VQSR_SNP_TRANCHES ;
    APPLY_VQSR ;
    ANNOTATE_FRAGMENTS
} from '../../../modules/local/vqsr'

// a resource VCF and its index: .tbi for bgzipped, .idx (Tribble) for plain text
def withIndex(path) {
    [file(path), file(path.endsWith('.gz') ? "${path}.tbi" : "${path}.idx")]
}

def stagedResources(vqsr, names) {
    names.collectMany { name -> withIndex(vqsr.resources[name]) }
}

// the shard names in a vcf export's manifest, in genomic order
def manifestNames(manifest) {
    manifest.readLines().collect { line -> line.trim() }.findAll { line -> line }
}

// a process output's files as a list: a glob matching one file gives a lone path, which is itself iterable
def asList(files) {
    files instanceof java.nio.file.Path ? [files] : files as List
}

// a vcf export, in chunks: [[VDS name, chunk index, shard VCFs, their indexes]], in genomic order
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

    COMPOSE_SITES_ONLY(densified.map { plan, _sites_only, separate_header, _whole -> [plan, separate_header] })
    def composed = COMPOSE_SITES_ONLY.out.map { plan -> [plan, file(plan.composed)] }

    TRAIN_VQSR_INDEL_MODEL(composed, indel_refs, vqsr)
    TRAIN_VQSR_SNP_MODEL(composed, snp_refs, vqsr)

    def chunks = densified.flatMap { plan, sites_only, _separate_header, _whole ->
        chunked(plan.name, plan.sites_only, manifestNames(sites_only), vqsr.fragments_per_job)
    }

    // the whole fragments, in chunks of their own size. Whole fragments and sites-only shards are both exported from
    // the dense MT's partitions, and listed in partition order, so the nth whole fragment holds the same records as the
    // nth sites-only shard. Their names aren't expected to match: they're paired by position alone, so only their
    // counts are checked
    def whole_chunks = densified.flatMap { plan, sites_only, _separate_header, whole ->
        def names = manifestNames(whole)
        def shard_count = manifestNames(sites_only).size()
        if (names.size() != shard_count) {
            error("${plan.name}: ${names.size()} whole fragments but ${shard_count} sites-only shards, so they can't be paired")
        }
        chunked(plan.name, plan.whole_fragments, names, vqsr.annotate_fragments_per_job)
    }

    // [VDS name, its number of sites-only shards], to size the regrouped VQSR VCFs' chunks
    def shard_counts = densified.map { plan, sites_only, _separate_header, _whole -> [plan.name, manifestNames(sites_only).size()] }

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

    // each VQSR VCF alone, by its sites-only shard's position in the export: [VDS name, position, VCF, index].
    // APPLY_VQSR names each after its input shard, but a glob output's order is its own, so they're matched by name
    def vqsr_shards = chunks
        .join(APPLY_VQSR.out.shards, by: [0, 1])
        .flatMap { key, index, shards, _tbis, vcfs, tbis ->
            def vcf_by_name = asList(vcfs).collectEntries { vcf -> [vcf.name, vcf] }
            def tbi_by_name = asList(tbis).collectEntries { tbi -> [tbi.name, tbi] }
            def names = shards.collect { shard -> "${shard.name.replaceFirst(/\.bgz$/, '')}.vqsr.vcf.gz".toString() }
            def missing = names.findAll { name -> !vcf_by_name.containsKey(name) || !tbi_by_name.containsKey("${name}.tbi".toString()) }
            if (missing) {
                error("${key} chunk ${index}: APPLY_VQSR wrote no ${missing}")
            }
            (0..<names.size()).collect { offset ->
                def name = names[offset]
                [key, index * vqsr.fragments_per_job + offset, vcf_by_name[name], tbi_by_name["${name}.tbi".toString()]]
            }
        }

    // the VQSR VCFs regrouped into the whole fragments' chunks, in position order: [VDS name, chunk index, VCFs,
    // indexes]. Each group's size is known, so it's emitted once its last VCF is, not once every APPLY_VQSR task is
    def annotate_size = vqsr.annotate_fragments_per_job
    def vqsr_annotate_chunks = vqsr_shards
        .combine(shard_counts, by: 0)
        .map { key, position, vcf, tbi, total ->
            def chunk = position.intdiv(annotate_size)
            [groupKey([key, chunk], Math.min(annotate_size, total - chunk * annotate_size)), position, vcf, tbi]
        }
        .groupTuple()
        .map { group, positions, vcfs, tbis ->
            def target = group.getGroupTarget()
            def order = (0..<positions.size()).toSorted { i -> positions[i] }
            [target[0], target[1], order.collect { i -> vcfs[i] }, order.collect { i -> tbis[i] }]
        }

    // each chunk of whole fragments meets the VQSR VCFs of the same partitions, both in partition order
    ANNOTATE_FRAGMENTS(whole_chunks.join(vqsr_annotate_chunks, by: [0, 1]))

    // every chunk writes the same header: chunk 0's is used
    def headers = ANNOTATE_FRAGMENTS.out.header
        .filter { _key, index, _header -> index == 0 }
        .map { key, _index, text -> [key, text] }

    emit:
    // [VDS name, chunk index, its annotated whole fragments, each a VCF with the same header]: chunks in no order
    fragments = ANNOTATE_FRAGMENTS.out.fragments
    // [VDS name, the annotated whole fragments' header, as text]
    header = headers
}
