/*
 * Combining for one dataset and sequencing type: Metamist cohorts (or a sample sheet) -> gVCFs -> one Hail VDS,
 * built on an ephemeral Dataproc cluster on top of the newest earlier VDS, and registered in Metamist.
 */

include { BUILD_SAMPLESHEET } from './modules/local/build_samplesheet'
include { COMBINE } from './modules/local/combine'

// per-run records: the sample sheet, the cluster and job, and what the combine added and removed
def runPath(name) {
    "combiner_samplesheets/${params.arGuid}/${name}"
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

// CLI params arrive as strings, and Groovy treats the string 'false' as true
def flag(param) {
    param.toString().toBoolean()
}

def intervals() {
    params.sequencing_type == 'exome' ? params.exome_intervals : params.genome_intervals
}

def validateParams() {
    def required = [
        'dataset', 'sequencing_type', 'arGuid', 'cpg_container',
        'dataproc_service_account', 'dataproc_staging_bucket', 'dataproc_autoscaling_policy',
    ]
    def missing = required.findAll { name -> !params[name] }
    if (missing) {
        error("Missing required params: ${missing.collect { name -> '--' + name }.join(', ')}")
    }
    if (!(params.sequencing_type in ['genome', 'exome'])) {
        error("--sequencing_type must be genome or exome, got '${params.sequencing_type}'")
    }
    if (!(params.accessLevel in ['full', 'standard', 'test'])) {
        error("--accessLevel must be full, standard or test, got '${params.accessLevel}'")
    }
    if ([params.cohorts, params.samplesheet].count { input -> input } != 1) {
        error('Provide exactly one of --cohorts or --samplesheet')
    }
    // a supplied sheet is an alternative input, and comes with an alternative output: its own VDS path
    if (params.samplesheet as boolean != params.output_vds as boolean) {
        error('--samplesheet and --output_vds go together: a sample-sheet run writes to the VDS path it is given')
    }
    if (params.output_vds && !params.output_vds.toString().endsWith('.vds')) {
        error("--output_vds must end in .vds, got '${params.output_vds}'")
    }
    ['force', 'register', 'preemptible_secondary_workers'].each { name ->
        if (!(params[name].toString() in ['true', 'false'])) {
            error("--${name} must be true or false, got '${params[name]}'")
        }
    }
    if (flag(params.force) && params.base_vds) {
        error('--force makes a fresh build, so it cannot be combined with --base_vds')
    }
    if (!intervals()) {
        error("No intervals for ${params.sequencing_type}: set --${params.sequencing_type}_intervals")
    }
}

// a sample sheet as a list of row maps, failing on a missing column, an empty value or a repeated SG
def readSheet(sheet) {
    def columns = ['sg_id', 'gvcf', 'dataset']
    def lines = sheet.readLines().findAll { line -> line.trim() }
    def header = lines ? lines[0].split('\t').toList() : []
    if (!header.containsAll(columns)) {
        fail("${sheet.name} needs the columns ${columns}, found ${header}")
    }
    def rows = lines.drop(1).collect { line -> [header, line.split('\t', -1).toList()].transpose().collectEntries() }
    def incomplete = rows.findAll { row -> columns.any { column -> !row[column] } }
    if (incomplete) {
        fail("${incomplete.size()} sample sheet row(s) have empty values, e.g. ${incomplete[0]}")
    }
    def repeated = rows.countBy { row -> row.sg_id }.findAll { _id, n -> n > 1 }.keySet()
    if (repeated) {
        fail("SGs appear more than once in the sample sheet: ${repeated.toSorted()}")
    }
    if (!rows) {
        fail('The sample sheet has no SGs')
    }
    return rows
}

/*
 * A cheap check on the VDS existence, all other checks are done at SampleSheet construction.
 */
def planCombine(sheet, output, base_vds, refresh, force) {
    def rows = readSheet(sheet)
    def sg_ids = rows.collect { row -> row.sg_id }.toSorted()

    // check the SUCCESS marker inside the VDS for existence, not the directory
    def exists = file("${output}/variant_data/_SUCCESS").exists()

    return [
        name: output.tokenize('/').last().replaceAll(/\.vds$/, ''),
        output: output,
        base_vds: base_vds,
        refresh: refresh.toSorted(),
        force: force,
        n_sgs: sg_ids.size(),
        action: exists ? 'done' : 'combine',
    ]
}

// error() inside a channel operator reaches the console only as "Unexpected error", so log the reason first
def fail(message) {
    log.error(message)
    error(message)
}

workflow {
    main:
    validateParams()
    def force = flag(params.force)

    // tolerate a list in config, or a ','-delimited String on CLI
    def refresh = params.refresh ? idList(params.refresh).unique() : []
    def cohorts = params.cohorts ? idList(params.cohorts).unique() : []

    // only cohort runs register: a sample-sheet run has nothing to link the VDS to
    def register = flag(params.register) && cohorts

    // [sheet, output VDS, base VDS or null]
    def built_sheet = channel.empty()
    def inputs = channel.empty()
    if (params.samplesheet) {
        // a test or one-off run, with its own output path. Metamist isn't consulted, so there's no base VDS
        // unless --base_vds names one, and no VDS counts as registered
        def output = uriOf(file(params.output_vds))
        inputs = channel.value([file(params.samplesheet, checkIfExists: true), output, null])
    }
    else {
        // the sheet is named for the VDS it makes (ADR 0004): <hash of SG and refresh IDs>_<SG count>.tsv
        def built = BUILD_SAMPLESHEET(metamistProject(), params.sequencing_type, cohorts, refresh)
        built_sheet = built.map { sheet, _base -> sheet }
        inputs = built.map { sheet, base ->
            def output = uriOf(workflow.outputDir.resolve("vds/${params.sequencing_type}/${sheet.baseName}.vds"))
            [sheet, output, base.text.trim() ?: null]
        }
    }

    def planned = inputs
        .map { sheet, output, base ->
            [sheet, planCombine(sheet, output, base, refresh, force)]
        }
        .branch { _sheet, plan ->
            done: plan.action == 'done'
            combine: true
        }

    // still gotta remove registered from here
    planned.done.subscribe { _sheet, plan ->
        log.info("${plan.output} already exists and is registered: these ${plan.n_sgs} SGs are already combined")
    }

    def bucket = "gs://cpg-${params.dataset}-${params.accessLevel == 'test' ? 'test' : 'main'}"
    def cluster = [
        ar_guid: params.arGuid,
        dataset: params.dataset,
        project: params.dataproc_project ?: params.dataset,
        region: params.dataproc_region,
        service_account: params.dataproc_service_account,
        staging_bucket: params.dataproc_staging_bucket,
        temp_bucket: "${bucket}-tmp",
        autoscaling_policy: params.dataproc_autoscaling_policy,
        master_type: params.dataproc_master_type,
        worker_type: params.dataproc_worker_type,
        num_workers: params.dataproc_num_workers,
        num_secondary_workers: params.dataproc_num_secondary_workers,
        preemptible_secondary_workers: flag(params.preemptible_secondary_workers),
        boot_disk_size_gb: params.dataproc_boot_disk_size_gb,
        num_local_ssds: params.dataproc_num_local_ssds,
        max_age: params.dataproc_max_age,
        max_idle: params.dataproc_max_idle,
        intervals: intervals(),
        branch_factor: params.branch_factor,
        target_records: params.target_records,
        gvcf_batch_size: params.gvcf_batch_size,
    ]
    def registration = register
        ? [project: metamistProject(), cohorts: cohorts, sequencing_type: params.sequencing_type, ar_guid: params.arGuid]
        : [:]

    COMBINE(planned.combine, cluster, registration)

    publish:
    samplesheet = built_sheet
    cluster = COMBINE.out.cluster
    report = COMBINE.out.report
}

// paths are relative to outputDir, which defaults to the dataset bucket (nextflow.config)
output {
    samplesheet {
        path { sheet -> sheet >> runPath('samplesheet.tsv') }
    }
    cluster {
        path { json -> json >> runPath('output.json') }
    }
    report {
        path { json -> json >> runPath('combine_report.json') }
    }
}
