/*
 * Densified sites-only exports -> one VQSR-filtered sites-only VCF per VDS, in the temp bucket.
 *
 * densified: [plan from planVqsr() in main.nf, header-per-shard manifest, separate-header manifest,
 *            genotypes-inclusive header-per-shard manifest]. The last isn't used here, only passed through to the emit
 * vqsr: settings from main.nf. fragments_per_job, the filter levels, and resources: {name: VCF path}
 *
 * The models are trained on the whole callset, joined by COMPOSE_SITES_ONLY. The SNP tranches are trained and
 * applied on chunks of header-per-shard shards. Everything is keyed by the VDS name, as a run can't tell how many
 * VDSs it will densify until it's planned them.
 */

include { COMPOSE_SITES_ONLY } from '../../../modules/local/compose'
include {
    TRAIN_VQSR_INDEL_MODEL ;
    TRAIN_VQSR_SNP_MODEL ;
    TRAIN_VQSR_SNP_TRANCHES ;
    GATHER_VQSR_SNP_TRANCHES ;
    APPLY_VQSR_SNP ;
    CONCAT_VCFS ;
    APPLY_VQSR_INDEL
} from '../../../modules/local/vqsr'

// a resource VCF and its index: .tbi for bgzipped, .idx (Tribble) for plain text
def withIndex(path) {
    [file(path), file(path.endsWith('.gz') ? "${path}.tbi" : "${path}.idx")]
}

def stagedResources(vqsr, names) {
    names.collectMany { name -> withIndex(vqsr.resources[name]) }
}

workflow VQSR {
    take:
    densified
    vqsr

    main:
    def snp_refs = stagedResources(vqsr, ['hapmap', 'omni', 'one_thousand_genomes', 'dbsnp'])
    def indel_refs = stagedResources(vqsr, ['mills', 'axiom_poly', 'dbsnp'])

    COMPOSE_SITES_ONLY(densified.map { plan, _sites_only, separate_header, _whole_fragments -> [plan, separate_header] })
    def composed = COMPOSE_SITES_ONLY.out.map { plan -> [plan, file(plan.composed)] }

    TRAIN_VQSR_INDEL_MODEL(composed, indel_refs, vqsr)
    TRAIN_VQSR_SNP_MODEL(composed, snp_refs, vqsr)

    // [VDS name, chunk index, shard VCFs, their indexes], the shards in manifest (genomic) order
    def chunks = densified.flatMap { plan, sites_only, _separate_header, _whole_fragments ->
        def names = sites_only.readLines().collect { line -> line.trim() }.findAll { line -> line }
        def batches = names.collate(vqsr.fragments_per_job)
        (0..<batches.size()).collect { index ->
            def shards = batches[index].collect { name -> "${plan.sites_only}/${name}" }
            [plan.name, index, shards.collect { shard -> file(shard) }, shards.collect { shard -> file("${shard}.tbi") }]
        }
    }

    TRAIN_VQSR_SNP_TRANCHES(chunks.combine(TRAIN_VQSR_SNP_MODEL.out, by: 0), snp_refs, vqsr)

    // every shard's tranches, gathered once all chunks are trained
    def shard_tranches = TRAIN_VQSR_SNP_TRANCHES.out.tranches
        .groupTuple()
        .map { key, tranches -> [key, tranches.flatten().toSorted { tranche -> tranche.name }] }
    GATHER_VQSR_SNP_TRANCHES(shard_tranches)

    // each chunk meets its own shards' recalibrations, and the gathered tranches
    def to_apply = chunks
        .join(TRAIN_VQSR_SNP_TRANCHES.out.recal, by: [0, 1])
        .combine(GATHER_VQSR_SNP_TRANCHES.out, by: 0)
    APPLY_VQSR_SNP(to_apply, vqsr)

    // the chunks back in genomic order, for a plain concat
    def snp_chunks = APPLY_VQSR_SNP.out
        .map { key, index, vcf, tbi -> [key, [index, vcf, tbi]] }
        .groupTuple()
        .map { key, parts ->
            def ordered = parts.toSorted { part -> part[0] }
            [key, ordered.collect { part -> part[1] }, ordered.collect { part -> part[2] }]
        }
    CONCAT_VCFS(snp_chunks)

    def to_indel = TRAIN_VQSR_INDEL_MODEL.out
        .map { plan, recal, recal_idx, tranches -> [plan.name, plan, recal, recal_idx, tranches] }
        .join(CONCAT_VCFS.out)
        .map { _key, plan, recal, recal_idx, tranches, vcf, tbi -> [plan, vcf, tbi, recal, recal_idx, tranches] }
    APPLY_VQSR_INDEL(to_indel, vqsr)

    // the genotypes-inclusive manifest rejoins its VDS's result; its directory is plan.whole_fragments
    def whole_manifests = densified.map { plan, _sites_only, _separate_header, whole_fragments -> [plan.name, whole_fragments] }
    def result = APPLY_VQSR_INDEL.out
        .map { plan, vcf, tbi -> [plan.name, plan, vcf, tbi] }
        .join(whole_manifests)
        .map { _key, plan, vcf, tbi, whole_fragments -> [plan, vcf, tbi, whole_fragments] }

    emit:
    // [plan, VQSR sites-only VCF, its index, genotypes-inclusive header-per-shard manifest]
    vcf = result
}
