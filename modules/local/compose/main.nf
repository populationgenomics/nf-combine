/*
 * The separate-header sites-only export -> one VCF, joined in GCS by `gcloud storage objects compose`
 * (bin/write_compose_script.py), so nothing is downloaded. plan comes from planVqsr() in main.nf, and manifest is
 * DENSIFY's copy of the export's shard manifest.
 *
 * Port of cpg-flow-seqr-loader's ConcatenateVcfsWithGcloud. The joined VCF is in the temp bucket, beside the
 * fragments: a compose can't cross buckets. A finished compose writes its output in one step, so if the output
 * exists it's complete, and is reused.
 */

process COMPOSE_SITES_ONLY {
    tag "${plan.name}"
    container params.cpg_container

    input:
    tuple val(plan), path(manifest)

    output:
    val plan

    script:
    """
    # something in the Batch runtime sets CLOUDSDK_PYTHON=/usr/bin/python3, absent in this image; let gcloud find its own
    unset CLOUDSDK_PYTHON

    if gcloud storage objects describe ${plan.composed} > /dev/null 2>&1; then
        echo "${plan.composed} exists, reusing it"
        exit 0
    fi

    # sizes, to skip the fragments of empty partitions
    gcloud storage ls -l ${plan.separate_header}/ > fragment_sizes.txt

    write_compose_script.py \\
        --manifest ${manifest} \\
        --vcf_dir ${plan.separate_header} \\
        --sizes fragment_sizes.txt \\
        --output ${plan.composed} \\
        --temp_dir ${plan.compose_temp} \\
        --script compose.sh

    bash compose.sh
    """

    stub:
    """
    echo "would compose ${plan.separate_header} into ${plan.composed}"
    """
}

/*
 * Fragments already in GCS -> one object, joined by `gcloud storage objects compose`. fragments are gs:// URIs in
 * the order they're joined. They're named, not staged, so nothing is downloaded. A compose can't cross buckets, so
 * the joining happens in compose_dir, which must be in the fragments' bucket. If the plan has a final_output, the
 * joined object is then copied there. A copy writes its output in one step, so if final_output exists it's complete,
 * and is reused.
 *
 * Every fragment is composed, empty or not: an empty object adds nothing. Used to join the VQSR VCF (VQSR).
 */
process COMPOSE_FRAGMENTS {
    tag "${plan.name}"
    container params.cpg_container

    input:
    tuple val(plan), val(fragments), val(compose_dir)

    output:
    // [plan, where the joined object ends up]
    tuple val(plan), val(plan.final_output ?: "${compose_dir}/${plan.name}.vcf.bgz")

    script:
    def reuse = plan.final_output
        ? """
        if gcloud storage objects describe ${plan.final_output} > /dev/null 2>&1; then
            echo "${plan.final_output} exists, reusing it"
            exit 0
        fi
        """
        : ''
    def final_output = plan.final_output ? "--final_output ${plan.final_output}" : ''
    """
    # something in the Batch runtime sets CLOUDSDK_PYTHON=/usr/bin/python3, absent in this image; let gcloud find its own
    unset CLOUDSDK_PYTHON
    ${reuse}
    # printf is a builtin, so the line length isn't limited
    printf '%s\\n' ${fragments.join(' ')} > fragments.txt

    write_compose_script.py \\
        --manifest fragments.txt \\
        --output ${compose_dir}/${plan.name}.vcf.bgz \\
        --temp_dir ${compose_dir}/intermediates \\
        ${final_output} \\
        --script compose.sh

    bash compose.sh
    """

    stub:
    """
    echo "would compose ${fragments.size()} fragments in ${compose_dir}, then copy to ${plan.final_output}"
    """
}
