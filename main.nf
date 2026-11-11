/*
 * Combining for one dataset and sequencing type: Metamist cohorts -> gVCFs -> one Hail VDS,
 * built on an ephemeral Dataproc cluster on top of the newest earlier VDS, and registered in Metamist.
 * Unless its VQSR VCF, or every dataset's VCF of it, already exists, the VDS is then densified, on a second cluster,
 * into sites-only VCF shards, and those are filtered with VQSR.
 * VQSR's results are copied onto the whole (genotypes-inclusive) VCF, which is joined into one VCF of every sample:
 * the VQSR VCF, kept in the temp bucket. Each dataset's VCF is subset from it, and registered in the dataset's
 * Metamist project. Datasets in seqr_vcf also get a copy of theirs with the samples named by participant, not SG.
 */

include { BUILD_SAMPLESHEET } from './modules/local/build_samplesheet'
include { COMBINE } from './modules/local/combine'
include { DENSIFY } from './modules/local/densify'
include { VQSR } from './subworkflows/local/vqsr'
include { SUBSET_DATASET_VCF ; REHEADER_VCF ; REGISTER_DATASET_VCF } from './modules/local/dataset_vcf'

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

/*
 * The VQSR VCF of a VDS: every sample, with VQSR's filters and annotations, and what each dataset's VCF is subset
 * from. In the temp bucket (<output_root>/tmp for local runs: Nextflow checks it exists, and can't read the bucket
 * offline), named for the VDS, so a later run finds it. Once the bucket's lifecycle rules remove it, a missing dataset
 * VCF means densifying again.
 */
def vqsrVcf(name) {
    def root = params.output_root ? "${uriOf(file(params.output_root))}/tmp/nf-combine" : tempRoot()
    "${root}/vqsr/${name}.vcf.bgz"
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
 * Also checks whether the VDS's VQSR VCF, and every dataset's VCF of it, exist: densifying and VQSR only make the
 * VQSR VCF, and it's only for the dataset VCFs, so if either is done, the VDS isn't densified.
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
        vqsr_vcf: vqsrVcf(sheet.baseName),
        // a compose or copy writes its output in one step, so if the VQSR VCF exists it's complete
        vqsr_vcf_done: file(vqsrVcf(sheet.baseName)).exists(),
        dataset_vcfs_done: planDatasetVcfs(sheet).every { _vds, plan -> plan.done },
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
 * Where VQSR on a densified VDS composes the sites-only VCF its models are trained on: in the temp bucket, named for
 * the VDS. It must share a bucket with the separate-header shards it's composed from. VQSR's product is the VQSR VCF
 * (vqsrVcf()), joined from the annotated whole fragments.
 */
def planVqsr(densify_plan) {
    def name = densify_plan.name
    def dir = "${tempRoot()}/vqsr/${name}"
    return densify_plan + [
        composed: "${dir}/${name}_sites_only.vcf.bgz",
        compose_temp: "${dir}/compose_intermediates",
        vqsr_vcf: vqsrVcf(name),
    ]
}

/*
 * One VCF per dataset with SGs in the sheet: the VQSR VCF, cut down to the dataset's SGs and the sites where they
 * have a non-reference allele. Written to the dataset's bucket (<output_root>/<dataset> for local runs), named for
 * the VDS. Its index is uploaded first, and an upload writes its object in one step, so if the VCF exists, an earlier
 * run finished it. Returns [[VDS name, plan]].
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
        def output = "${root}/joint_vcf/${params.sequencing_type}/${vds}.vcf.bgz"
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
            def seqr_output = "${root}/joint_vcf/${params.sequencing_type}/${vds}_seqr.vcf.bgz"
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

    // a VDS is densified, whether made by this run or already there, unless its VQSR VCF, or every dataset's VCF of
    // it, exists
    def vds_plans = planned.done.map { _sheet, plan -> plan }
        .mix(COMBINE.out.plan)
        .branch { plan ->
            dataset_vcfs_done: plan.dataset_vcfs_done
            vqsr_vcf_done: plan.vqsr_vcf_done
            todo: true
        }

    vds_plans.dataset_vcfs_done.subscribe { plan ->
        log.info("Every dataset VCF of ${plan.output} already exists: it isn't densified")
    }
    vds_plans.vqsr_vcf_done.subscribe { plan ->
        log.info("${plan.vqsr_vcf} already exists: ${plan.output} isn't densified, and the dataset VCFs are subset from it")
    }

    DENSIFY(vds_plans.todo.map { plan -> planDensify(plan) }, hail_lib, cluster)

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

    VQSR(densified, vqsr)

    // [VDS name, its VQSR VCF]: made by this run, or already there. Without --vqsr, no dataset VCFs are made
    def vqsr_vcfs = vds_plans.vqsr_vcf_done
        .filter { _plan -> params.vqsr }
        .map { plan -> [plan.name, plan.vqsr_vcf] }
        .mix(VQSR.out.vcf)

    def dataset_plans = built
        .flatMap { sheet, _base -> planDatasetVcfs(sheet) }
        .branch { _vds, plan ->
            done: plan.done
            todo: true
        }

    dataset_plans.done.subscribe { _vds, plan ->
        log.info("${plan.final_output} already exists: the ${plan.dataset} VCF is already done")
    }

    // each missing dataset VCF, from its VDS's VQSR VCF
    SUBSET_DATASET_VCF(dataset_plans.todo.combine(vqsr_vcfs, by: 0).map { _vds, plan, vcf -> [plan, vcf] })

    // only the VCFs this run made: an existing one was registered by the run that made it
    REGISTER_DATASET_VCF(SUBSET_DATASET_VCF.out.filter { _plan, _vcf -> params.register }, registration)

    // seqr VCFs: from each seqr_vcf dataset's VCF, made by this run or already there
    def seqr_plans = dataset_plans.done
        .map { _vds, plan -> [plan, plan.final_output] }
        .mix(SUBSET_DATASET_VCF.out)
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
