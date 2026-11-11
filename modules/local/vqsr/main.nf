/*
 * GATK VQSR on the sites-only VCFs, ported from cpg-flow-seqr-loader's TrainVqsr*, GatherTrainedVqsrSnpTranches,
 * RunSnpVqsrOnFragments and RunIndelVqsr stages. The models are trained on the composed VCF (COMPOSE_SITES_ONLY).
 * The SNP tranches are trained shard by shard, on the header-per-shard export, with many shards per task. SNP and
 * indel VQSR are then applied to each shard, and each shard's VQSR VCF is copied onto the whole fragment from the same
 * partition (ANNOTATE_FRAGMENTS), ready to be joined into the VQSR VCF.
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

/*
 * bash: strip_eof FILE cuts the BGZF EOF block off the end of FILE, after checking it's there. bcftools ends every
 * bgzipped file with one, and htsjdk (GATK) reads one in the middle of a file as the end of the file, so it's cut
 * from everything that's composed but the final eof.bgz
 */
def stripEof() {
    '''
    strip_eof() {
        if [[ $(tail -c 28 "$1" | od -An -tx1 | tr -d ' \\n') != 1f8b08040000000000ff0600424302001b0003000000000000000000 ]]; then
            echo "$1 doesn't end in a BGZF EOF block" >&2
            exit 1
        fi
        truncate -s $(( $(wc -c < "$1") - 28 )) "$1"
    }
    '''
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
 * One chunk of shards: SNP then indel VQSR applied to each. Each shard's VQSR VCF is kept: ANNOTATE_FRAGMENTS copies
 * it onto the whole fragment from the same partition. The indel model is genome-wide, and ApplyVQSR looks its
 * recalibrations up by locus, so applying it shard by shard gives the same result as applying it to the joined VCF.
 */
process APPLY_VQSR {
    tag "${key} chunk ${chunk}"
    container params.gatk_container

    input:
    tuple val(key), val(chunk), path(vcfs), path(tbis), path(snp_recals), path(snp_recal_idxs), path(snp_tranches), path(indel_recal), path(indel_recal_idx), path(indel_tranches)
    val vqsr

    output:
    tuple val(key), val(chunk), path('*.vqsr.vcf.gz'), path('*.vqsr.vcf.gz.tbi'), emit: shards

    script:
    """
    for vcf in ${vcfs}; do
        ${linkShard()}
        gatk --java-options "${javaOptions(task)}" ApplyVQSR \\
            --tmp-dir . \\
            -V \$name.vcf.bgz \\
            -O \$name.snp.vcf.gz \\
            --recal-file \$name.recal \\
            --tranches-file ${snp_tranches} \\
            --truth-sensitivity-filter-level ${vqsr.snp_filter_level} \\
            --use-allele-specific-annotations \\
            -mode SNP
        gatk --java-options "${javaOptions(task)}" ApplyVQSR \\
            --tmp-dir . \\
            -V \$name.snp.vcf.gz \\
            -O \$name.vqsr.vcf.gz \\
            --recal-file ${indel_recal} \\
            --tranches-file ${indel_tranches} \\
            --truth-sensitivity-filter-level ${vqsr.indel_filter_level} \\
            --use-allele-specific-annotations \\
            -mode INDEL
        rm \$name.snp.vcf.gz \$name.snp.vcf.gz.tbi
    done
    """

    stub:
    """
    for vcf in ${vcfs}; do
        touch \${vcf%.bgz}.vqsr.vcf.gz \${vcf%.bgz}.vqsr.vcf.gz.tbi
    done
    """
}

/*
 * One chunk of whole (genotypes-inclusive) fragments: each gets the FILTER and VQSR INFO fields of the VQSR VCF
 * from the same partition, and is written ready to compose into the VQSR VCF (COMPOSE_FRAGMENTS): headless, and
 * without a BGZF EOF block. The chunk's header is written alone, bgzipped and without an EOF block, and so is the EOF
 * block, to go first and last in the compose.
 *
 * Each fragment is annotated against its own shard's VQSR VCF, not one joined across shards: without a region,
 * bcftools reads the annotation file from the start of the fragment's chromosome, so every fragment would read most
 * of a chromosome. Fragments and VQSR VCFs come in partition order, and are paired by position: their names aren't
 * expected to match. Each pair must hold the same records: their record counts are compared before annotating.
 *
 * Every chunk's header gets VQSR's header lines from -h, so every chunk writes the same header. An empty partition's
 * fragment is an empty file: composing it adds nothing.
 */
process ANNOTATE_FRAGMENTS {
    tag "${key} chunk ${chunk}"
    container params.cpg_container

    input:
    tuple val(key), val(chunk), path(fragments), path(fragment_tbis), path(vqsr_vcfs), path(vqsr_tbis)

    output:
    tuple val(key), val(chunk), path('*.annotated.vcf.bgz'), emit: fragments
    tuple val(key), val(chunk), path('header.vcf.bgz'), path('eof.bgz'), emit: header

    script:
    """
    ${stripEof()}

    # INFO IDs in a VCF's header, sorted for comm
    info_ids() {
        bcftools view -h "\$1" | sed -n 's/^##INFO=<ID=\\([^,]*\\),.*/\\1/p' | LC_ALL=C sort
    }

    # both in partition order, as staged: the nth fragment pairs with the nth VQSR VCF
    fragments=(${fragments})
    vqsr_vcfs=(${vqsr_vcfs})
    if (( \${#fragments[@]} != \${#vqsr_vcfs[@]} )); then
        echo "\${#fragments[@]} fragments but \${#vqsr_vcfs[@]} VQSR VCFs: the exports' partitions don't match" >&2
        exit 1
    fi
    first=\${vqsr_vcfs[0]}

    # VQSR's INFO fields are the ones its VCF has and the fragments don't. The rest are left as the fragments have them
    vqsr_fields=\$(LC_ALL=C comm -13 <(info_ids \${fragments[0]}) <(info_ids \$first))
    if ! grep -qx AS_VQSLOD <<< "\$vqsr_fields"; then
        echo "\$first has no AS_VQSLOD the fragments lack: is it a VQSR VCF? New fields: \$vqsr_fields" >&2
        exit 1
    fi
    columns=FILTER\$(printf ',INFO/%s' \$vqsr_fields)
    echo "Copying \$columns"

    # VQSR's FILTER lines and its new INFO lines. Built from header lines alone, as any shard may be empty
    bcftools view -h \$first \\
        | grep -E "^##FILTER=|^##INFO=<ID=(\$(paste -sd '|' - <<< "\$vqsr_fields"))," > vqsr_header_lines.txt
    # the VQSR VCF's header comes from these lines alone
    if ! grep -q '^##INFO=<ID=AS_VQSLOD,' vqsr_header_lines.txt; then
        echo "No AS_VQSLOD header line was taken from \$first" >&2
        exit 1
    fi
    bcftools view --no-version -h \${fragments[0]} \\
        | bcftools annotate --no-version -h vqsr_header_lines.txt -Oz -o header.vcf.bgz
    tail -c 28 header.vcf.bgz > eof.bgz
    strip_eof header.vcf.bgz

    for i in \${!fragments[@]}; do
        vcf=\${fragments[i]}
        vqsr_vcf=\${vqsr_vcfs[i]}
        name=\${vcf%.bgz}
        # counted by reading the records: neither Hail's nor GATK's tabix indexes hold count metadata, so
        # `bcftools index -n` gives 0 for both, and every fragment would be taken as empty
        records=\$(bcftools view -H \$vcf | wc -l)
        if [[ \$records != \$(bcftools view -H \$vqsr_vcf | wc -l) ]]; then
            echo "\$vcf and \$vqsr_vcf hold different numbers of records: the exports' partitions don't match" >&2
            exit 1
        fi
        # nothing to annotate, and bcftools (1.21 at least) segfaults annotating an empty VCF from another
        if (( records == 0 )); then
            : > \$name.annotated.vcf.bgz
            continue
        fi
        # exact: a record takes VQSR's fields only from the record with all of its alleles. Written headless, as the
        # header is composed in once, first
        bcftools annotate --no-version --pair-logic exact -a \$vqsr_vcf -c \$columns -h vqsr_header_lines.txt \\
            --threads ${task.cpus} -Ou \$vcf \\
            | bcftools view --no-version -H --threads ${task.cpus} -Oz -o \$name.annotated.vcf.bgz
        strip_eof \$name.annotated.vcf.bgz
    done
    """

    stub:
    """
    touch header.vcf.bgz eof.bgz
    for vcf in ${fragments}; do
        touch \${vcf%.bgz}.annotated.vcf.bgz
    done
    """
}
