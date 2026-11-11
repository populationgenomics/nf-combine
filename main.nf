/*
 * Combining for one dataset and sequencing type: Metamist cohorts -> gVCFs -> one Hail VDS,
 * built on an ephemeral Dataproc cluster on top of the newest earlier VDS, and registered in Metamist.
 * The VDS is then densified, on a second cluster, into sites-only VCF shards, and those are filtered with VQSR.
 * VQSR's results are copied onto the whole (genotypes-inclusive) VCF, which is then split into one VCF per dataset.
 * Each dataset's VCF is registered in its Metamist project. Datasets in seqr_vcf also get a copy of theirs with the
 * samples named by participant, not SG.
 */

include { BUILD_SAMPLESHEET } from './modules/local/build_samplesheet'
include { COMBINE } from './modules/local/combine'
include { DENSIFY } from './modules/local/densify'
include { VQSR } from './subworkflows/local/vqsr'
include { DATASET_VCFS } from './subworkflows/local/dataset_vcfs'
include { REHEADER_VCF ; REGISTER_DATASET_VCF } from './modules/local/dataset_vcf'

/*
 * Types only: the defaults are in nextflow.config, which conf/base.config reads before this script is parsed.
 * Nextflow converts command-line values to these types (`--force false` is a Boolean), and fails the run if a
 * param without `?` has no value.
 */
params {
    dataset: String
    sequencing_type: String
    accessLevel: String
    arGuid: String
    cohorts: String  // or a list, from a config file

    force: Boolean
    base_vds: String?
    refresh: String? // or a list, from a config file
    register: Boolean
    output_root: String?
    seqr_vcf: String? // or a list, from a config file

    // only the run's sequencing type's are needed: validateParams() checks them
    exome_combiner_intervals: String?
    exome_densify_intervals: String?
    genome_combiner_intervals: String?
    genome_densify_intervals: String?

    target_records: Integer
    branch_factor: Integer
    gvcf_batch_size: Integer

    dataproc_project: String
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

    vqsr: Boolean
    vqsr_fragments_per_job: Integer
    annotate_fragments_per_job: Integer
    snp_filter_level: Float
    indel_filter_level: Float
    vqsr_hapmap: String
    vqsr_omni: String
    vqsr_one_thousand_genomes: String
    vqsr_dbsnp: String
    vqsr_mills: String
    vqsr_axiom_poly: String

    cpg_container: String
    gatk_container: String
}

/*
 * Output destinations, worked out here and handed to the processes that publish them: publishDir
 * (conf/base.config) only applies them, so the VDS path and the run records can't drift apart.
 */

// a dataset's bucket, at the run's access level
def bucketOf(dataset) {
    "gs://cpg-${dataset}-${params.accessLevel == 'test' ? 'test' : 'main'}"
}

// the workflow-running dataset's bucket
def bucket() {
    bucketOf(params.dataset)
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

// a dataset's Metamist project, at the run's access level
def metamistProjectOf(dataset) {
    params.accessLevel == 'test' ? "${dataset}-test" : dataset
}

def metamistProject() {
    metamistProjectOf(params.dataset)
}

// gs://bucket/path for cloud paths (Path.toString() drops the scheme), a plain path locally. No trailing slash: a
// cloud directory's toUriString() ends in one
def uriOf(path) {
    (path.scheme == 'file' ? path.toString() : path.toUriString()).replaceAll('/+$', '')
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

// the datasets whose VCF is also written with participant IDs as sample names
def seqrDatasets() {
    params.seqr_vcf ? idList(params.seqr_vcf).unique() : []
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
 * The sheet's participant IDs become the sample names of the seqr VCFs, so in each seqr_vcf dataset every SG needs
 * one, and no two SGs can share one. Checked before anything is combined. A seqr_vcf dataset with no SGs in the
 * sheet gets no VCF, so it's warned about.
 */
def checkSeqrParticipants(sheet) {
    def by_dataset = sheet.splitCsv(header: true, sep: '\t').groupBy { row -> row.dataset }
    seqrDatasets().findAll { dataset -> !by_dataset.containsKey(dataset) }.each { dataset ->
        log.warn("--seqr_vcf names ${dataset}, which has no SGs in the sample sheet: it gets no VCF")
    }
    def errors = []
    by_dataset.subMap(seqrDatasets()).each { dataset, rows ->
        def missing = rows.findAll { row -> !row.participant }.collect { row -> row.sg_id }
        if (missing) {
            errors << "${dataset}: no participant ID for ${missing}"
        }
        rows.findAll { row -> row.participant }.groupBy { row -> row.participant }.each { participant, shared ->
            if (shared.size() > 1) {
                errors << "${dataset}: participant ${participant} has several SGs: ${shared.collect { row -> row.sg_id }}"
            }
        }
    }
    if (errors) {
        error("Can't name the seqr VCF samples by participant:\n${errors.join('\n')}")
    }
}

/*
 * A cheap check on the VDS existence, all other checks are done at SampleSheet construction.
 */
def planCombine(sheet, discovered_base, refresh, force) {
    // the sheet is named for the VDS it makes: <hash of SG and refresh IDs>_<SG count>.tsv
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
        // fragments of a genotypes-inclusive VCF representation
        whole_fragments: "${dir}/${name}_full.vcf.bgz",
        run_dir: combine_plan.run_dir,
    ]
}

/*
 * Where VQSR on a densified VDS writes, in the temp bucket beside the densified outputs and named for the VDS.
 * The composed VCF must share a bucket with the separate-header shards it's composed from.
 */
def planVqsr(densify_plan) {
    def name = densify_plan.name
    def dir = "${tempRoot()}/vqsr/${name}"
    def vqsr_name = "${name}_vqsr.vcf.gz"
    return densify_plan + [
        composed: "${dir}/${name}_sites_only.vcf.bgz",
        compose_temp: "${dir}/compose_intermediates",
        vqsr_dir: dir,
        vqsr_name: vqsr_name,
        vqsr_vcf: "${dir}/${vqsr_name}",
        // the index is published with the VCF, so if it exists, an earlier run finished VQSR on this VDS
        done: file("${dir}/${vqsr_name}.tbi").exists(),
    ]
}

/*
 * One VCF per dataset with SGs in the sheet: the VQSR-annotated whole VCF, cut down to the dataset's SGs and the
 * sites where they have a non-reference allele. Written to the dataset's bucket (<output_root>/<dataset> for local
 * runs), named for the VDS. If it exists, an earlier run finished it. Returns [[VDS name, plan]].
 * metamist_project is the dataset's own project, where its VCF is registered.
 *
 * A seqr_vcf dataset's plan also has seqr_output, beside the VCF: the same VCF, with each SG's sample renamed to
 * its participant ID (checkSeqrParticipants()). Its index is written last, so if that exists, the VCF is done.
 */
def planDatasetVcfs(sheet) {
    def vds = sheet.baseName
    def by_dataset = sheet.splitCsv(header: true, sep: '\t').groupBy { row -> row.dataset }
    by_dataset.collect { dataset, rows ->
        def root = params.output_root ? "${uriOf(file(params.output_root))}/${dataset}" : bucketOf(dataset)
        def output = "${root}/vcf/${params.sequencing_type}/${vds}.vcf.bgz"
        def plan = [
            name: "${vds}_${dataset}",
            dataset: dataset,
            metamist_project: metamistProjectOf(dataset),
            // sorted, as the sample order is the VCF's
            sgs: rows.collect { row -> row.sg_id }.toSorted(),
            final_output: output,
            done: file(output).exists(),
        ]
        if (dataset in seqrDatasets()) {
            def seqr_output = "${root}/vcf/${params.sequencing_type}/${vds}_seqr.vcf.bgz"
            plan += [
                seqr_output: seqr_output,
                participants: rows.collectEntries { row -> [row.sg_id, row.participant] },
                seqr_done: file("${seqr_output}.tbi").exists(),
            ]
        }
        [vds, plan]
    }
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
            checkSeqrParticipants(sheet)
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
        // a dataset's GCP project ID can differ from its name (seqr is seqr-nnnnn)
        dataset: params.dataset,
        project: params.dataproc_project,
        region: params.dataproc_region,
        service_account: "seqera-${params.accessLevel}@${params.dataproc_project}.iam.gserviceaccount.com",
        temp_bucket: tempBucket(),
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

    def vqsr = [
        fragments_per_job: params.vqsr_fragments_per_job,
        annotate_fragments_per_job: params.annotate_fragments_per_job,
        snp_filter_level: params.snp_filter_level,
        indel_filter_level: params.indel_filter_level,
        resources: [
            hapmap: params.vqsr_hapmap,
            omni: params.vqsr_omni,
            one_thousand_genomes: params.vqsr_one_thousand_genomes,
            dbsnp: params.vqsr_dbsnp,
            mills: params.vqsr_mills,
            axiom_poly: params.vqsr_axiom_poly,
        ],
    ]
    def densified = DENSIFY.out.shards
        .filter { _plan, _sites_only, _separate_header, _whole_fragments -> params.vqsr }
        .map { plan, sites_only, separate_header, whole_fragments -> [planVqsr(plan), sites_only, separate_header, whole_fragments] }
        .branch { plan, _sites_only, _separate_header, _whole_fragments ->
            done: plan.done
            todo: true
        }

    densified.done.subscribe { plan, _sites_only, _separate_header, _whole_fragments ->
        log.info("${plan.vqsr_vcf} already exists: VQSR on this VDS is already done")
    }

    VQSR(densified.todo, vqsr)

    def dataset_plans = built
        .flatMap { sheet, _base -> planDatasetVcfs(sheet) }
        .branch { _vds, plan ->
            done: plan.done
            todo: true
        }

    dataset_plans.done.subscribe { _vds, plan ->
        log.info("${plan.final_output} already exists: the ${plan.dataset} VCF is already done")
    }

    DATASET_VCFS(VQSR.out.fragments, VQSR.out.header, dataset_plans.todo)

    // only the VCFs this run made: an existing one was registered by the run that made it
    REGISTER_DATASET_VCF(DATASET_VCFS.out.vcf.filter { _plan, _vcf -> params.register }, registration)

    // seqr VCFs: from each seqr_vcf dataset's VCF, made by this run or already there
    def seqr_plans = dataset_plans.done
        .map { _vds, plan -> [plan, plan.final_output] }
        .mix(DATASET_VCFS.out.vcf)
        .filter { plan, _vcf -> plan.seqr_output }
        .branch { plan, _vcf ->
            done: plan.seqr_done
            todo: true
        }

    seqr_plans.done.subscribe { plan, _vcf ->
        log.info("${plan.seqr_output} already exists: the ${plan.dataset} seqr VCF is already done")
    }

    REHEADER_VCF(seqr_plans.todo)
}
