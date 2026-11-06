/*
 * GATK VQSR on the sites-only VCFs, ported from cpg-flow-seqr-loader's TrainVqsr*, GatherTrainedVqsrSnpTranches,
 * RunSnpVqsrOnFragments and RunIndelVqsr stages. The models are trained on the composed VCF (COMPOSE_SITES_ONLY).
 * The SNP tranches are trained and applied shard by shard, on the header-per-shard export, with many shards per
 * task. Indel VQSR is applied last, to the joined SNP-recalibrated VCF.
 *
 * vqsr holds the settings (planned in main.nf): filter levels, and the resources, as {name: VCF path}. The
 * resource files and their indexes are staged as refs, and are found by their basenames.
 */

// VariantRecalibrator's arguments in either mode, apart from its input and outputs: as in cpg-flow-seqr-loader
def recalibratorArgs(mode, vqsr) {
    def snp = mode == 'SNP'
    // depth (AS_VarDP) is left out, as in the Broad pipelines: it's meaningless for exomes
    def features = ['AS_ReadPosRankSum', 'AS_MQRankSum', 'AS_QD', 'AS_FS', 'AS_SOR'] + (snp ? ['AS_MQ'] : [])
    def tranches = snp
        ? [100.0, 99.95, 99.9, 99.8, 99.6, 99.5, 99.4, 99.3, 99.0, 98.0, 97.0, 90.0]
        : [100.0, 99.95, 99.9, 99.5, 99.0, 97.0, 96.0, 95.0, 94.0, 93.5, 93.0, 92.0, 91.0, 90.0]
    // each training resource's role and prior
    def resources = snp
        ? [
            hapmap: 'known=false,training=true,truth=true,prior=15',
            omni: 'known=false,training=true,truth=true,prior=12',
            one_thousand_genomes: 'known=false,training=true,truth=false,prior=10',
            dbsnp: 'known=true,training=false,truth=false,prior=7',
        ]
        : [
            mills: 'known=false,training=true,truth=true,prior=12',
            axiom_poly: 'known=false,training=true,truth=false,prior=10',
            dbsnp: 'known=true,training=false,truth=false,prior=2',
        ]
    def args = ['--trust-all-polymorphic', "-mode ${mode}", '--use-allele-specific-annotations']
    args += tranches.collect { value -> "-tranche ${value}" }
    args += features.collect { feature -> "-an ${feature}" }
    args += resources.collect { name, role -> "-resource:${name},${role} ${vqsr.resources[name].tokenize('/').last()}" }
    args.join(' \\\n        ')
}

// heap at 80% of the task's memory, leaving the rest for the JVM and GATK's native code
def javaOptions(task) {
    "-Xmx${(task.memory.toMega() * 4).intdiv(5)}m -XX:+UseParallelGC -XX:ParallelGCThreads=${task.cpus}"
}

// bash: links shard $vcf as $name.vcf.bgz. GATK picks a VCF codec by extension, and Hail names shards part-NNNNN.bgz
def linkShard() {
    'name=${vcf%.bgz}; ln -sf $vcf $name.vcf.bgz; ln -sf $vcf.tbi $name.vcf.bgz.tbi'
}

process TRAIN_VQSR_INDEL_MODEL {
    tag "${plan.name}"
    container params.gatk_container

    input:
    tuple val(plan), path(sites_only)
    path refs
    val vqsr

    output:
    tuple val(plan), path('indel.recal'), path('indel.recal.idx'), path('indel.tranches')

    script:
    """
    gatk IndexFeatureFile -I ${sites_only}

    # --max-gaussians 4: indels are fewer than SNPs. A small callset can still fail with 'No data found', if so lower it
    gatk --java-options "${javaOptions(task)}" VariantRecalibrator \\
        -V ${sites_only} \\
        -O indel.recal \\
        --tranches-file indel.tranches \\
        --max-gaussians 4 \\
        ${recalibratorArgs('INDEL', vqsr)}
    """

    stub:
    """
    touch indel.recal indel.recal.idx indel.tranches
    """
}

process TRAIN_VQSR_SNP_MODEL {
    tag "${plan.name}"
    container params.gatk_container

    input:
    tuple val(plan), path(sites_only)
    path refs
    val vqsr

    output:
    tuple val(plan.name), path('snp.model')

    script:
    """
    gatk IndexFeatureFile -I ${sites_only}

    # only the model is kept: the tranches are trained shard by shard, against it
    gatk --java-options "${javaOptions(task)}" VariantRecalibrator \\
        -V ${sites_only} \\
        -O snp.recal \\
        --tranches-file snp.tranches \\
        --output-model snp.model \\
        --sample-every-Nth-variant 10 \\
        --max-gaussians 6 \\
        ${recalibratorArgs('SNP', vqsr)}
    """

    stub:
    """
    touch snp.model
    """
}

/*
 * One chunk of header-per-shard shards, scored against the SNP model. A chunk shares one image pull and one set of
 * staged resources: per shard, those were nearly all of the cost. Each shard's outputs are named for it.
 */
process TRAIN_VQSR_SNP_TRANCHES {
    tag "${key} chunk ${chunk}"
    container params.gatk_container

    input:
    tuple val(key), val(chunk), path(vcfs), path(tbis), path(snp_model)
    path refs
    val vqsr

    output:
    tuple val(key), val(chunk), path('*.recal'), path('*.recal.idx'), emit: recal
    tuple val(key), path('*.tranches'), emit: tranches

    script:
    """
    for vcf in ${vcfs}; do
        ${linkShard()}
        gatk --java-options "${javaOptions(task)}" VariantRecalibrator \\
            --verbosity WARNING \\
            --QUIET \\
            -V \$name.vcf.bgz \\
            -O \$name.recal \\
            --tranches-file \$name.tranches \\
            --input-model ${snp_model} \\
            --output-tranches-for-scatter \\
            --max-gaussians 6 \\
            ${recalibratorArgs('SNP', vqsr)}
        # as in cpg-flow, in case the index isn't written in scatter mode
        touch \$name.recal.idx
    done
    """

    stub:
    """
    for vcf in ${vcfs}; do
        touch \${vcf%.bgz}.recal \${vcf%.bgz}.recal.idx \${vcf%.bgz}.tranches
    done
    """
}

process GATHER_VQSR_SNP_TRANCHES {
    tag "${key}"
    container params.gatk_container

    input:
    tuple val(key), path(tranches)

    output:
    tuple val(key), path('snp.tranches')

    script:
    def inputs = tranches.collect { tranche -> "--input ${tranche}" }.join(' ')
    """
    gatk --java-options "${javaOptions(task)}" GatherTranches --mode SNP ${inputs} --output snp.tranches
    """

    stub:
    """
    touch snp.tranches
    """
}

/*
 * One chunk of shards: SNP VQSR applied to each, then the chunk gathered into one VCF. The shards arrive in
 * genomic order, so the chunk's VCF is sorted.
 */
process APPLY_VQSR_SNP {
    tag "${key} chunk ${chunk}"
    container params.gatk_container

    input:
    tuple val(key), val(chunk), path(vcfs), path(tbis), path(recals), path(recal_idxs), path(tranches)
    val vqsr

    output:
    tuple val(key), val(chunk), path("chunk_${chunk}.vcf.gz"), path("chunk_${chunk}.vcf.gz.tbi")

    script:
    """
    gather=()
    for vcf in ${vcfs}; do
        ${linkShard()}
        gatk --java-options "${javaOptions(task)}" ApplyVQSR \\
            -V \$name.vcf.bgz \\
            -O \$name.vqsr.vcf.gz \\
            --recal-file \$name.recal \\
            --tranches-file ${tranches} \\
            --truth-sensitivity-filter-level ${vqsr.snp_filter_level} \\
            --use-allele-specific-annotations \\
            -mode SNP
        gather+=(-I \$name.vqsr.vcf.gz)
    done

    gatk --java-options "${javaOptions(task)}" GatherVcfs "\${gather[@]}" -O chunk_${chunk}.vcf.gz
    gatk IndexFeatureFile -I chunk_${chunk}.vcf.gz
    """

    stub:
    """
    touch chunk_${chunk}.vcf.gz chunk_${chunk}.vcf.gz.tbi
    """
}

// sorted, non-overlapping VCFs, in order -> one VCF
process CONCAT_VCFS {
    tag "${key}"
    container params.bcftools_container

    input:
    tuple val(key), path(vcfs), path(tbis)

    output:
    tuple val(key), path('snp_vqsr.vcf.gz'), path('snp_vqsr.vcf.gz.tbi')

    script:
    """
    bcftools concat --threads ${task.cpus} ${vcfs} -Oz -o snp_vqsr.vcf.gz -W=tbi
    """

    stub:
    """
    touch snp_vqsr.vcf.gz snp_vqsr.vcf.gz.tbi
    """
}

// published to plan.vqsr_dir (conf/base.config), for the annotation stages to read
process APPLY_VQSR_INDEL {
    tag "${plan.name}"
    container params.gatk_container

    input:
    tuple val(plan), path(snp_vcf), path(snp_tbi), path(recal), path(recal_idx), path(tranches)
    val vqsr

    output:
    tuple val(plan), path(plan.vqsr_name), path("${plan.vqsr_name}.tbi")

    script:
    """
    gatk --java-options "${javaOptions(task)}" ApplyVQSR \\
        --tmp-dir . \\
        -V ${snp_vcf} \\
        -O ${plan.vqsr_name} \\
        --recal-file ${recal} \\
        --tranches-file ${tranches} \\
        --truth-sensitivity-filter-level ${vqsr.indel_filter_level} \\
        --use-allele-specific-annotations \\
        -mode INDEL
    """

    stub:
    """
    touch ${plan.vqsr_name} ${plan.vqsr_name}.tbi
    """
}
