/*
 * One dataset's share of the VQSR-annotated whole VCF: its SGs, and only the sites where at least one of them has a
 * non-reference allele. The header and the fragments are subset in separate tasks, so a dataset's header is made
 * without staging any fragment, and both are written ready to compose (COMPOSE_FRAGMENTS): the header, then the
 * headless fragments, then the BGZF EOF block.
 *
 * plan comes from planDatasetVcfs() in main.nf: its sgs are the dataset's SGs, the samples kept.
 */

// bcftools view's options for a dataset's subset, shared by its header and its fragments, so they match. --min-ac 1
// drops the rows where none of the samples has a non-reference allele, counted from their genotypes. The export is
// split (densify.py), one alt per row, so a row is kept or dropped whole: no allele needs trimming. --no-update
// keeps every INFO field, AC and AN included, as computed across the whole callset
def subsetArgs() {
    '--no-version --no-update --samples-file samples.txt --min-ac 1'
}

// bash: the dataset's SGs, one per line. The order is the VCF's sample order, so the same everywhere. printf is a
// builtin, so the line length isn't limited
def writeSamples(plan) {
    "printf '%s\\n' ${plan.sgs.join(' ')} > samples.txt"
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

// the annotated fragments' header (text) -> the dataset's header, bgzipped without an EOF block, and the EOF block
process SUBSET_HEADER {
    tag "${plan.name}"
    container params.bcftools_container

    input:
    tuple val(plan), path(header)

    output:
    tuple val(plan.name), path('header.vcf.bgz'), path('eof.bgz')

    script:
    """
    ${stripEof()}
    ${writeSamples(plan)}

    bcftools view ${subsetArgs()} -h ${header} -Oz -o header.vcf.bgz
    tail -c 28 header.vcf.bgz > eof.bgz
    strip_eof header.vcf.bgz
    """

    stub:
    """
    touch header.vcf.bgz eof.bgz
    """
}

// one chunk of annotated fragments -> the dataset's subset of each, headless, without an EOF block. A fragment with
// no sites left is empty
process SUBSET_FRAGMENTS {
    tag "${plan.name} chunk ${chunk}"
    container params.bcftools_container

    input:
    tuple val(plan), val(chunk), path(fragments)

    output:
    tuple val(plan.name), val(chunk), path('*.subset.vcf.bgz')

    script:
    """
    ${stripEof()}
    ${writeSamples(plan)}

    for vcf in ${fragments}; do
        name=\${vcf%.annotated.vcf.bgz}
        bcftools view ${subsetArgs()} -H --threads ${task.cpus} -Oz -o \$name.subset.vcf.bgz \$vcf
        strip_eof \$name.subset.vcf.bgz
    done
    """

    stub:
    """
    for vcf in ${fragments}; do
        touch \${vcf%.annotated.vcf.bgz}.subset.vcf.bgz
    done
    """
}
