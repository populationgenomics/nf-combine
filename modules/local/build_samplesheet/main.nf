// cohorts -> one combining sample sheet, named for the VDS it makes, plus the base VDS to build on;
// see bin/build_samplesheet.py
process BUILD_SAMPLESHEET {
    tag "${project}"

    input:
    val project
    val sequencing_type
    val cohorts
    val refresh

    output:
    tuple path('*.tsv'), path('base_vds.txt')

    script:
    def refresh_ids = refresh ? "--refresh ${refresh.join(' ')}" : ''
    """
    build_samplesheet.py \\
        --project ${project} \\
        --sequencing-type ${sequencing_type} \\
        --cohorts ${cohorts.join(' ')} \\
        ${refresh_ids}
    """

    // task.ext.fixtures, set from config in the stub tests, is a directory holding both outputs
    stub:
    def fixtures = task.ext.fixtures
    if (fixtures) {
        """
        cp ${fixtures}/*.tsv ${fixtures}/base_vds.txt .
        """
    }
    else {
        """
        printf 'sg_id\\tgvcf\\tdataset\\n' > stub_0.tsv
        touch base_vds.txt
        """
    }
}
