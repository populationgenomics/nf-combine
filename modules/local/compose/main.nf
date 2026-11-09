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
    container params.gcloud_container

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
