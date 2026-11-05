/*
 * One combine on an ephemeral Dataproc cluster, then registration of the VDS (ADR 0003). plan comes from
 * planCombine() in main.nf; cluster holds the Dataproc and Hail settings. hail_lib is the nf_combine_hail package,
 * shipped to the cluster for combiner.py to import. registration is empty for runs that don't register, else
 * [project, cohorts, sequencing_type, ar_guid].
 */
include { clusterArgs } from '../dataproc'

process COMBINE {
    tag "${plan.name}"

    input:
    tuple path(samplesheet), val(plan)
    path hail_lib
    val cluster
    val registration

    output:
    val plan, emit: plan
    path 'output.json', emit: cluster
    path 'combine_report.json', emit: report

    script:
    def base = plan.base_vds ? "--base_vds ${plan.base_vds}" : ''
    def refresh = plan.refresh ? "--refresh ${plan.refresh.join(' ')}" : ''
    def register = registration ? registerCommand(plan, registration, '') : ''
    """
    # bin/ is on PATH, and the Hail script is uploaded to the cluster from there. A failed lookup inside the
    # argument list would pass silently, so it's looked up on its own line
    HAIL_SCRIPT=\$(which combiner.py)
    submit_hail.py \\
        --cluster-name combine-${cluster.ar_guid} \\
        ${clusterArgs(cluster)} \\
        --script \$HAIL_SCRIPT \\
        --pyfiles ${hail_lib} \\
        --files ${samplesheet} \\
        --complete ${plan.output}/reference_data/_SUCCESS ${plan.output}/variant_data/_SUCCESS \\
        --fetch ${plan.report} combine_report.json \\
        --output output.json \\
        --record vds=${plan.output} base_vds=${plan.base_vds ?: ''} \\
        -- \\
        --samplesheet ${samplesheet} \\
        --out ${plan.output} \\
        --temp_dir ${plan.temp_dir} \\
        --report ${plan.report} \\
        --intervals ${cluster.combiner_intervals} \\
        --branch ${cluster.branch_factor} \\
        --target_records ${cluster.target_records} \\
        --gvcf_batch ${cluster.gvcf_batch_size} ${base} ${refresh} ${plan.force ? '--force' : ''}

    ${register}
    """

    stub:
    def register = registration ? registerCommand(plan, registration, '--dry') : ''
    """
    printf '{"vds": "${plan.output}", "base_vds": "${plan.base_vds ?: ''}", "cluster": null}\\n' > output.json
    printf '{"vds": "${plan.output}", "added": [], "removed": [], "refreshed": []}\\n' > combine_report.json
    echo ${register}
    """
}

// the VDS is linked to the cohorts; base_vds and refresh are recorded only when used
def registerCommand(plan, registration, dry) {
    def meta = ["sequencing_type=${registration.sequencing_type}", "ar-guid=${registration.ar_guid}"]
    if (plan.base_vds) {
        meta << "base_vds=${plan.base_vds}"
    }
    if (plan.refresh) {
        meta << "refresh=${plan.refresh.join(',')}"
    }
    "python3 -m cpg_utils.metamist_registration ${dry} --project ${registration.project} --output ${plan.output} " +
        "--type vds --cohorts ${registration.cohorts.join(' ')} --meta ${meta.join(' ')}"
}
