/*
 * The VQSR VCF -> one dataset's VCF: its SGs, and only the sites where at least one of them has a non-reference
 * allele, indexed, at plan.final_output. plan comes from planDatasetVcfs() in main.nf: its sgs are the dataset's SGs,
 * the samples kept, sorted. bcftools writes them in the order listed, so that's the VCF's sample order.
 *
 * The VQSR VCF is copied into a scratch directory on the task's disk (not the work directory, which would keep it),
 * and everything there is deleted when the task ends, whether or not it succeeds. The index is uploaded first, then
 * the VCF: an upload writes its object in one step, so if the VCF exists, it's complete and indexed
 * (planDatasetVcfs()). A dataset VCF written before they were indexed also counts as done.
 *
 * --min-ac 1 drops the rows where none of the samples has a non-reference allele, counted from their genotypes. The
 * export is split (densify.py), one alt per row, so a row is kept or dropped whole: no allele needs trimming.
 * --no-update keeps every INFO field, AC and AN included, as computed across the whole callset.
 */
process SUBSET_DATASET_VCF {
    tag "${plan.name}"
    container params.cpg_container
    // the VQSR VCF, the dataset's share of it (no bigger), and 20 GB to spare
    disk { 1.GB * ((2 * file(vcf).size()).intdiv(1000000000) + 21) }

    input:
    // vcf: the VQSR VCF, as a gs:// URI. Not staged: it's copied into scratch, then deleted
    tuple val(plan), val(vcf)

    output:
    tuple val(plan), val(plan.final_output)

    script:
    """
    # something in the Batch runtime sets CLOUDSDK_PYTHON=/usr/bin/python3, absent in this image; let gcloud find its own
    unset CLOUDSDK_PYTHON

    scratch=\$(mktemp -d)
    trap 'rm -rf "\$scratch"' EXIT

    # printf is a builtin, so the line length isn't limited
    printf '%s\\n' ${plan.sgs.join(' ')} > samples.txt

    gcloud storage cp ${vcf} \$scratch/input.vcf.bgz
    bcftools view --no-version --no-update --samples-file samples.txt --min-ac 1 \\
        --threads ${task.cpus} -Oz -o \$scratch/dataset.vcf.bgz \$scratch/input.vcf.bgz
    rm \$scratch/input.vcf.bgz

    # exactly the dataset's SGs, in order
    bcftools query -l \$scratch/dataset.vcf.bgz > subset_samples.txt
    if ! diff subset_samples.txt samples.txt; then
        echo "${plan.name}: the subset VCF's samples aren't the dataset's SGs" >&2
        exit 1
    fi

    bcftools index --tbi --threads ${task.cpus} \$scratch/dataset.vcf.bgz
    gcloud storage cp \$scratch/dataset.vcf.bgz.tbi ${plan.final_output}.tbi
    gcloud storage cp \$scratch/dataset.vcf.bgz ${plan.final_output}
    """

    stub:
    """
    echo "would subset ${vcf} to ${plan.sgs.size()} SGs, to ${plan.final_output}"
    """
}

/*
 * A dataset's VCF -> the same VCF with each SG's sample renamed to its participant's primary external ID (from the
 * sample sheet), indexed, at plan.seqr_output. Only for the datasets in params.seqr_vcf.
 *
 * The VCF is copied into a scratch directory on the task's disk (not the work directory, which would keep it), and
 * everything there is deleted when the task ends, whether or not it succeeds. bcftools reheader only rewrites the
 * header and copies the rest of the file block by block, but it writes a new file, so the disk holds two copies.
 * The index is uploaded last: if it exists, the VCF is complete (planDatasetVcfs() in main.nf).
 */
process REHEADER_VCF {
    tag "${plan.name}"
    container params.cpg_container
    // the VCF twice (before and after reheadering), and 20 GB to spare
    disk { 1.GB * ((2 * file(vcf).size()).intdiv(1000000000) + 21) }

    input:
    // vcf: the dataset's VCF, as a gs:// URI. Not staged: it's copied into scratch, then deleted
    tuple val(plan), val(vcf)

    output:
    tuple val(plan), val(plan.seqr_output)

    script:
    // SG -> participant pairs, and the participants in the VCF's sample order (plan.sgs), to check the result
    def pairs = plan.sgs.collect { sg -> "'${sg}' '${plan.participants[sg]}'" }.join(' ')
    def expected = plan.sgs.collect { sg -> "'${plan.participants[sg]}'" }.join(' ')
    """
    # something in the Batch runtime sets CLOUDSDK_PYTHON=/usr/bin/python3, absent in this image; let gcloud find its own
    unset CLOUDSDK_PYTHON

    scratch=\$(mktemp -d)
    trap 'rm -rf "\$scratch"' EXIT

    printf '%s\\t%s\\n' ${pairs} > participants.txt
    printf '%s\\n' ${expected} > expected_samples.txt

    gcloud storage cp ${vcf} \$scratch/input.vcf.bgz
    bcftools reheader --samples participants.txt --output \$scratch/seqr.vcf.bgz \$scratch/input.vcf.bgz
    rm \$scratch/input.vcf.bgz

    # every sample renamed, in the same order
    bcftools query -l \$scratch/seqr.vcf.bgz > samples.txt
    if ! diff samples.txt expected_samples.txt; then
        echo "${plan.name}: the reheadered VCF's samples aren't the expected participants" >&2
        exit 1
    fi

    bcftools index --tbi --threads ${task.cpus} \$scratch/seqr.vcf.bgz
    gcloud storage cp \$scratch/seqr.vcf.bgz ${plan.seqr_output}
    gcloud storage cp \$scratch/seqr.vcf.bgz.tbi ${plan.seqr_output}.tbi
    """

    stub:
    """
    echo "would reheader ${vcf} with participant IDs, to ${plan.seqr_output}"
    """
}

/*
 * Registers a dataset's VCF, made by this run, in the dataset's own Metamist project, as a dataset_vcf analysis linked
 * to every SG in it (plan.sgs). registration holds the run's sequencing_type and ar_guid. A retry only repeats the
 * registration: the VCF is already written.
 */
process REGISTER_DATASET_VCF {
    tag "${plan.name}"
    container params.cpg_container

    input:
    tuple val(plan), val(vcf)
    val registration

    output:
    tuple val(plan), val(vcf)

    script:
    registerDatasetVcf(plan, vcf, registration, '')

    stub:
    """
    echo ${registerDatasetVcf(plan, vcf, registration, '--dry')}
    """
}

def registerDatasetVcf(plan, vcf, registration, dry) {
    def meta = ["sequencing_type=${registration.sequencing_type}", "ar-guid=${registration.ar_guid}"]
    "python3 -m cpg_utils.metamist_registration ${dry} --project ${plan.metamist_project} --output ${vcf} " +
        "--type dataset_vcf --sgs ${plan.sgs.join(' ')} --meta ${meta.join(' ')}"
}
