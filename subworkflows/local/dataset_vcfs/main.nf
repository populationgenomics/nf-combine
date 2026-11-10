/*
 * VQSR-annotated whole fragments -> one VCF per dataset: its SGs, and the sites where at least one of them has a
 * non-reference allele.
 *
 * fragments: [VDS name, chunk index, annotated whole fragments], from VQSR
 * header: [VDS name, the annotated fragments' header, as text], from VQSR
 * datasets: [VDS name, dataset plan from planDatasetVcfs() in main.nf]
 *
 * Every chunk of fragments is subset for every dataset as soon as it's annotated, and each dataset's header is
 * subset on its own, from the header alone. A dataset's header, its headless fragments in genomic order, and the
 * BGZF EOF block are then joined in GCS (COMPOSE_FRAGMENTS), beside the dataset's header: every fragment is a task
 * output, so all are in the work directory's bucket. The joined VCF is then copied to the plan's final_output.
 */

include { SUBSET_HEADER ; SUBSET_FRAGMENTS } from '../../../modules/local/dataset_vcf'
include { COMPOSE_FRAGMENTS } from '../../../modules/local/compose'

// gs://bucket/path for cloud paths (Path.toString() drops the scheme), a plain path locally
def uriOf(path) {
    path.scheme == 'file' ? path.toString() : path.toUriString()
}

workflow DATASET_VCFS {
    take:
    fragments
    header
    datasets

    main:
    // [dataset VCF name, plan], to put the plan back on outputs keyed by name
    def plans = datasets.map { _vds, plan -> [plan.name, plan] }

    SUBSET_HEADER(header.combine(datasets, by: 0).map { _vds, text, plan -> [plan, text] })
    SUBSET_FRAGMENTS(fragments.combine(datasets, by: 0).map { _vds, index, vcfs, plan -> [plan, index, vcfs] })

    // the fragments in genomic order: chunks by index, and within one, by partition name (part-NNNNN, zero-padded).
    // A chunk of one fragment emits a lone path, hence the flatten. Then the header first and the EOF block last
    def to_compose = SUBSET_FRAGMENTS.out
        .map { name, index, vcfs -> [name, [index, [vcfs].flatten().toSorted { vcf -> vcf.name }]] }
        .groupTuple()
        .join(SUBSET_HEADER.out)
        .join(plans)
        .map { _name, parts, dataset_header, eof, plan ->
            def ordered = parts.toSorted { part -> part[0] }.collectMany { part -> part[1] }
            def joined = ([dataset_header] + ordered + [eof]).collect { vcf -> uriOf(vcf) }
            [plan, joined, "${uriOf(dataset_header.parent)}/compose"]
        }
    COMPOSE_FRAGMENTS(to_compose)

    emit:
    // [dataset plan, the dataset's VCF]
    vcf = COMPOSE_FRAGMENTS.out
}
