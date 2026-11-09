/*
 * One VDS -> a dense MT and two sites-only VCF exports, on an ephemeral Dataproc cluster (bin/densify.py). plan
 * comes from planDensify() in main.nf; cluster holds the Dataproc and Hail settings; hail_lib is the
 * nf_combine_hail package, shipped to the cluster for densify.py to import.
 *
 * Every output is in the temp bucket and named for the VDS, so a retry, or a later run on the same VDS, reuses
 * whatever is complete. Each export's shard manifest is copied into the task, for later processes to scatter on.
 */
include { clusterArgs } from '../dataproc'

process DENSIFY {
    tag "${plan.name}"

    input:
    val plan
    path hail_lib
    val cluster

    output:
    tuple val(plan), path('sites_only.manifest'), path('separate_header.manifest'), path('whole_fragments.manifest'), emit: shards
    path 'densify_output.json', emit: record

    script:
    """
    # a failed lookup inside the argument list would pass silently, so it's looked up on its own line
    HAIL_SCRIPT=\$(which densify.py)

    # the manifest files are fetched here so the next steps can scatter across them
    submit_hail.py \\
        --cluster-name densify-${cluster.ar_guid} \\
        ${clusterArgs(cluster)} \\
        --script \$HAIL_SCRIPT \\
        --pyfiles ${hail_lib} \\
        --complete ${plan.mt}/_SUCCESS ${plan.sites_only}/shard-manifest.txt ${plan.separate_header}/shard-manifest.txt ${plan.whole_fragments}/shard-manifest.txt \\
        --fetch ${plan.sites_only}/shard-manifest.txt sites_only.manifest \\
        --fetch ${plan.separate_header}/shard-manifest.txt separate_header.manifest \\
        --fetch ${plan.whole_fragments}/shard-manifest.txt whole_fragments.manifest \\
        --output densify_output.json \\
        --record vds=${plan.vds} mt=${plan.mt} sites_only=${plan.sites_only} separate_header=${plan.separate_header} whole_fragments=${plan.whole_fragments} \\
        -- \\
        --vds ${plan.vds} \\
        --intervals ${cluster.densify_intervals} \\
        --checkpoint ${plan.checkpoint} \\
        --out ${plan.mt} \\
        --sites_only ${plan.sites_only} \\
        --separate_header ${plan.separate_header} \\
        --whole_fragments ${plan.whole_fragments}
    """

    stub:
    """
    printf 'part-00000.bgz\\npart-00001.bgz\\n' > sites_only.manifest
    printf 'header.bgz\\npart-00000.bgz\\npart-00001.bgz\\n' > separate_header.manifest
    printf 'part-00000.bgz\\npart-00001.bgz\\n' > whole_fragments.manifest
    printf '{"vds": "${plan.vds}", "mt": "${plan.mt}", "cluster": null}\\n' > densify_output.json
    """
}
