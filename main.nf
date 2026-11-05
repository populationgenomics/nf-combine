/*
 * Combining for one dataset and sequencing type: Metamist cohorts -> gVCFs -> one Hail VDS,
 * built on an ephemeral Dataproc cluster on top of the newest earlier VDS, and registered in Metamist.
 * The VDS is then densified, on a second cluster, into sites-only VCF shards.
 */

include { BUILD_SAMPLESHEET } from './modules/local/build_samplesheet'
include { COMBINE } from './modules/local/combine'
include { DENSIFY } from './modules/local/densify'

/*
 * Types only: the defaults are in nextflow.config, which conf/base.config reads before this script is parsed.
 * Nextflow converts command-line values to these types (`--force false` is a Boolean), and fails the run if a
 * param without `?` has no value.
 */
params {
    dataset: String
    dataproc_project: String
    sequencing_type: String
    accessLevel: String
    arGuid: String
    cohorts: String  // or a list, from a config file

    force: Boolean
    base_vds: String?
    refresh: String? // or a list, from a config file
    register: Boolean
    output_root: String?

    // only the run's sequencing type's are needed: validateParams() checks them
    exome_combiner_intervals: String?
    exome_densify_intervals: String?
    genome_combiner_intervals: String?
    genome_densify_intervals: String?

    target_records: Integer
    branch_factor: Integer
    gvcf_batch_size: Integer

    dataproc_project: String?
    dataproc_region: String
    dataproc_autoscaling_policy: String
    dataproc_master_type: String
    dataproc_worker_type: String
    dataproc_num_workers: Integer
    dataproc_num_secondary_workers: Integer
    preemptible_secondary_workers: Boolean
    dataproc_boot_disk_size_gb: Integer
    dataproc_num_local_ssds: Integer
    dataproc_max_age: Integer
    dataproc_max_idle: Integer

    cpg_container: String
}

/*
 * Output destinations, worked out here and handed to the processes that publish them: publishDir
 * (conf/base.config) only applies them, so the VDS path and the run records can't drift apart.
 */

// the workflow-running dataset's bucket
def bucket() {
    "gs://cpg-${params.dataset}-${params.accessLevel == 'test' ? 'test' : 'main'}"
}

// the bucket, or output_root for local runs (absolute, whatever the launch dir)
def datasetRoot() {
    params.output_root ? uriOf(file(params.output_root)) : bucket()
}

// Dataproc's staging and temp bucket, and the root of the pipeline's own temp files: combiner plans and reports,
// and densified outputs. The bucket's lifecycle rules clean it up
def tempBucket() {
    "${bucket()}-tmp"
}

def tempRoot() {
    "${tempBucket()}/nf-combine"
}

// per-run records: the sample sheet, the cluster and job, and what the combine added and removed
def runDir() {
    "${datasetRoot()}/combiner_samplesheets/${params.arGuid}"
}

def metamistProject() {
    params.accessLevel == 'test' ? "${params.dataset}-test" : params.dataset
}

// gs://bucket/path for cloud paths (Path.toString() drops the scheme), a plain path locally
def uriOf(path) {
    path.scheme == 'file' ? path.toString() : path.toUriString()
}

// a list from config, or a comma/space separated string from the command line
def idList(param) {
    (param instanceof List) ? param : param.toString().tokenize(', ')
}

def combinerIntervals() {
    params.sequencing_type == 'exome' ? params.exome_combiner_intervals : params.genome_combiner_intervals
}

def densifyIntervals() {
    params.sequencing_type == 'exome' ? params.exome_densify_intervals : params.genome_densify_intervals
}

// what the params block can't check: allowed values, and params that are needed together or not at all
def validateParams() {
    if (!(params.sequencing_type in ['genome', 'exome'])) {
        error("--sequencing_type must be genome or exome, got '${params.sequencing_type}'")
    }
    if (!(params.accessLevel in ['full', 'standard', 'test'])) {
        error("--accessLevel must be full, standard or test, got '${params.accessLevel}'")
    }
    if (params.force && params.base_vds) {
        error('--force makes a fresh build, so it cannot be combined with --base_vds')
    }
    [combiner: combinerIntervals(), densify: densifyIntervals()].each { use, intervals ->
        def param = "--${params.sequencing_type}_${use}_intervals"
        if (!intervals) {
            error("No ${use} intervals set for ${params.sequencing_type}: set ${param}")
        }
        // read by Hail on the cluster, which can't see local files
        if (!intervals.startsWith('gs://')) {
            error("${param} must be a gs:// path, got '${intervals}'")
        }
    }
}

/*
 * A cheap check on the VDS existence, all other checks are done at SampleSheet construction.
 */
def planCombine(sheet, discovered_base, refresh, force) {
    // the sheet is named for the VDS it makes (ADR 0004): <hash of SG and refresh IDs>_<SG count>.tsv
    def output = "${datasetRoot()}/vds/${params.sequencing_type}/${sheet.baseName}.vds"

    // a supplied --base_vds is the user's to check; force always makes a fresh build
    def base_vds = force ? null : (params.base_vds ?: discovered_base)

    // check the SUCCESS marker inside the VDS for existence, not the directory
    def exists = file("${output}/variant_data/_SUCCESS").exists()

    return [
        name: sheet.baseName,
        output: output,
        base_vds: base_vds,
        refresh: refresh.toSorted(),
        force: force,
        n_sgs: sheet.readLines().count { line -> line.trim() } - 1,
        action: exists ? 'done' : 'combine',
        // keyed by the VDS name, so a retry finds the plan and report of its earlier attempt
        temp_dir: tempRoot(),
        report: "${tempRoot()}/combine_reports/${sheet.baseName}.json",
        // only read by publishDir (conf/base.config)
        run_dir: runDir(),
    ]
}

/*
 * Where densifying a VDS writes, all in the temp bucket and named for the VDS, so densifying the same VDS again
 * reuses what's complete. Nothing downstream reads the MT: it's kept only to resume from.
 */
def planDensify(combine_plan) {
    def name = combine_plan.name
    def dir = "${tempRoot()}/densify/${name}"
    return [
        name: name,
        vds: combine_plan.output,
        checkpoint: "${dir}/${name}_checkpoint.mt",
        mt: "${dir}/${name}.mt",
        // header per shard: each shard is a valid VCF, for processing shard by shard
        sites_only: "${dir}/${name}.vcf.bgz",
        // separate header: the shards concatenate, after the header, into one VCF
        separate_header: "${dir}/${name}_separate.vcf.bgz",
        run_dir: combine_plan.run_dir,
    ]
}

workflow {
    validateParams()

    // tolerate a list in config, or a ','-delimited String on CLI
    def refresh = params.refresh ? idList(params.refresh).unique() : []
    def cohorts = idList(params.cohorts).unique()

    // [sheet, base_vds.txt]
    def built = BUILD_SAMPLESHEET(metamistProject(), params.sequencing_type, cohorts, refresh, runDir())

    def planned = built
        .map { sheet, base ->
            [sheet, planCombine(sheet, base.text.trim() ?: null, refresh, params.force)]
        }
        .branch { _sheet, plan ->
            done: plan.action == 'done'
            combine: true
        }

    // subscribe is a non-blocking side-effect operator, so this only prints
    planned.done.subscribe { _sheet, plan ->
        log.info("${plan.output} already exists and is registered: these ${plan.n_sgs} SGs are already combined")
    }

    def cluster = [
        ar_guid: params.arGuid,
        // a dataset's GCP project ID can differ from its name (seqr is seqr-308602)
        project: params.dataproc_project ?: params.dataset,
        region: params.dataproc_region,
        service_account: "seqera-${params.accessLevel}@${params.dataproc_project}.iam.gserviceaccount.com",
        temp_dir: tempBucket(),
        autoscaling_policy: params.dataproc_autoscaling_policy,
        master_type: params.dataproc_master_type,
        worker_type: params.dataproc_worker_type,
        num_workers: params.dataproc_num_workers,
        num_secondary_workers: params.dataproc_num_secondary_workers,
        preemptible_secondary_workers: params.preemptible_secondary_workers,
        boot_disk_size_gb: params.dataproc_boot_disk_size_gb,
        num_local_ssds: params.dataproc_num_local_ssds,
        max_age: params.dataproc_max_age,
        max_idle: params.dataproc_max_idle,
        combiner_intervals: combinerIntervals(),
        densify_intervals: densifyIntervals(),
        branch_factor: params.branch_factor,
        target_records: params.target_records,
        gvcf_batch_size: params.gvcf_batch_size,
    ]

    def registration = params.register
        ? [project: metamistProject(), cohorts: cohorts, sequencing_type: params.sequencing_type, ar_guid: params.arGuid]
        : [:]

    // the Hail helpers the cluster's scripts import, shipped with each job
    def hail_lib = file("${projectDir}/hail_lib/nf_combine_hail")

    COMBINE(planned.combine, hail_lib, cluster, registration)

    // every VDS is densified, whether made by this run or already there
    def vds_plans = planned.done.map { _sheet, plan -> plan }.mix(COMBINE.out.plan)
    DENSIFY(vds_plans.map { plan -> planDensify(plan) }, hail_lib, cluster)
}
