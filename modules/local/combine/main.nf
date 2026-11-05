/*
 * One combine on an ephemeral Dataproc cluster, then registration of the VDS (ADR 0003). plan comes from
 * planCombine() in main.nf; cluster holds the Dataproc and combiner settings. registration is empty for runs
 * that don't register, else [project, cohorts, sequencing_type, ar_guid].
 */
process COMBINE {
    tag "${plan.name}"

    input:
    tuple path(samplesheet), val(plan)
    val cluster
    val registration

    output:
    path 'output.json', emit: cluster
    path 'combine_report.json', emit: report

    script:
    def base = plan.base_vds ? "--base-vds ${plan.base_vds}" : ''
    def refresh = plan.refresh ? "--refresh ${plan.refresh.join(' ')}" : ''
    def register = registration ? registerCommand(plan, registration, '') : ''
    """
    # bin/ is on PATH, and the Hail script is uploaded to the cluster from there. A failed lookup inside the
    # argument list would pass silently, so it's looked up on its own line
    HAIL_SCRIPT=\$(which combiner.py)
    submit_combiner.py \\
        --cluster-name combine-${cluster.ar_guid} \\
        --region ${cluster.region} \\
        --project ${cluster.project} \\
        --script \$HAIL_SCRIPT \\
        --samplesheet ${samplesheet} \\
        --staging-bucket ${cluster.staging_bucket} \\
        --temp-bucket ${cluster.temp_bucket} \\
        --vds ${plan.output} \\
        ${base} \\
        ${refresh} \\
        --intervals ${cluster.intervals} \\
        --autoscaling-policy ${cluster.autoscaling_policy} \\
        --master-type ${cluster.master_type} \\
        --worker-type ${cluster.worker_type} \\
        --num-workers ${cluster.num_workers} \\
        --num-secondary-workers ${cluster.num_secondary_workers} \\
        ${cluster.preemptible_secondary_workers ? '--preemptible-secondary-workers' : ''} \\
        --max-age ${cluster.max_age} \\
        --max-idle ${cluster.max_idle} \\
        --boot-disk-size ${cluster.boot_disk_size_gb} \\
        --num-local-ssds ${cluster.num_local_ssds} \\
        --service-account ${cluster.service_account} \\
        --branch-factor ${cluster.branch_factor} \\
        --target-records ${cluster.target_records} \\
        --gvcf-batch-size ${cluster.gvcf_batch_size} \\
        ${plan.force ? '--force' : ''} \\
        --labels ar-guid=${cluster.ar_guid} compute-category=dataproc dataset=${cluster.project}

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
