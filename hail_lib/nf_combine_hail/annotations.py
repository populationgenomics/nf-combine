"""Variant-level annotation expressions: strand bias, allele balance, low-quality and adj genotype flags.

Copied from cpg-flow-seqr-loader's hail_scripts/annotations.py (originally gnomad_methods), keeping only
what sparse_mt uses.
"""

import hail as hl


def to_phred(linear_expr: hl.expr.NumericExpression) -> hl.expr.Float64Expression:
    """Computes the phred-scaled value of a linear-scale input.

    Args:
        linear_expr: The linear-scale value.

    Returns:
        The phred-scaled value.
    """
    return -10 * hl.log10(linear_expr)


def get_lowqual_expr(
    alleles: hl.expr.ArrayExpression,
    qual_approx_expr: hl.expr.ArrayNumericExpression | hl.expr.NumericExpression,
    snv_phred_threshold: int = 30,
    snv_phred_het_prior: int = 30,  # 1/1000
    indel_phred_threshold: int = 30,
    indel_phred_het_prior: int = 39,  # 1/8,000
) -> hl.expr.BooleanExpression | hl.expr.ArrayExpression:
    """Computes the lowqual flag for split or unsplit alleles, from QUALapprox or AS_QUALapprox.

    Computed from QUALapprox, this differs from GATK's LowQual filter. GATK works at the site level, using the
    least stringent prior for mixed sites. Computed from AS_QUALapprox, this can be more stringent for some
    alleles at mixed sites.

    Args:
        alleles: The site's alleles.
        qual_approx_expr: QUALapprox, or the per-alt-allele AS_QUALapprox array.
        snv_phred_threshold: Phred-scaled SNV emission threshold, like GATK's.
        snv_phred_het_prior: Phred-scaled SNV heterozygosity prior (30 = 1/1000 bases, the GATK default).
        indel_phred_threshold: Phred-scaled indel emission threshold, like GATK's.
        indel_phred_het_prior: Phred-scaled indel heterozygosity prior.

    Returns:
        One flag per alt allele for AS_QUALapprox, else one flag for the site.
    """
    min_snv_qual = snv_phred_threshold + snv_phred_het_prior
    min_indel_qual = indel_phred_threshold + indel_phred_het_prior
    min_mixed_qual = max(min_snv_qual, min_indel_qual)

    if isinstance(qual_approx_expr, hl.expr.ArrayNumericExpression):
        return hl.range(1, hl.len(alleles)).map(
            lambda ai: hl.if_else(
                hl.is_snp(alleles[0], alleles[ai]),
                qual_approx_expr[ai - 1] < min_snv_qual,
                qual_approx_expr[ai - 1] < min_indel_qual,
            )
        )
    return (
        hl.case()
        .when(
            hl.range(1, hl.len(alleles)).all(lambda ai: hl.is_snp(alleles[0], alleles[ai])),
            qual_approx_expr < min_snv_qual,
        )
        .when(
            hl.range(1, hl.len(alleles)).all(lambda ai: hl.is_indel(alleles[0], alleles[ai])),
            qual_approx_expr < min_indel_qual,
        )
        .default(qual_approx_expr < min_mixed_qual)
    )


def get_adj_expr(
    gt_expr: hl.expr.CallExpression,
    gq_expr: hl.expr.Int32Expression | hl.expr.Int64Expression,
    dp_expr: hl.expr.Int32Expression | hl.expr.Int64Expression,
    ad_expr: hl.expr.ArrayNumericExpression,
    adj_gq: int = 20,
    adj_dp: int = 10,
    adj_ab: float = 0.2,
    haploid_adj_dp: int = 5,
) -> hl.expr.BooleanExpression:
    """Flags genotypes that pass the adj (adjusted) quality thresholds. The defaults are gnomAD's.

    Args:
        gt_expr: The genotype.
        gq_expr: Genotype quality.
        dp_expr: Read depth.
        ad_expr: Allele depths.
        adj_gq: Minimum GQ.
        adj_dp: Minimum DP for diploid genotypes.
        adj_ab: Minimum allele balance for each called allele of a het.
        haploid_adj_dp: Minimum DP for haploid genotypes.

    Returns:
        True where the genotype passes.
    """
    return (
        (gq_expr >= adj_gq)
        & hl.if_else(gt_expr.is_haploid(), dp_expr >= haploid_adj_dp, dp_expr >= adj_dp)
        & (
            hl.case()
            .when(~gt_expr.is_het(), True)
            .when(gt_expr.is_het_ref(), ad_expr[gt_expr[1]] / dp_expr >= adj_ab)
            .default((ad_expr[gt_expr[0]] / dp_expr >= adj_ab) & (ad_expr[gt_expr[1]] / dp_expr >= adj_ab))
        )
    )


def fs_from_sb(
    sb: hl.expr.ArrayNumericExpression | hl.expr.ArrayExpression,
    normalize: bool = True,
    min_cell_count: int = 200,
    min_count: int = 4,
    min_p_value: float = 1e-320,
) -> hl.expr.Int64Expression:
    """Computes FS (Fisher strand balance) from an SB (strand balance) table.

    FS is the phred-scaled p-value of a two-sided Fisher exact test on strand balance. With the defaults this
    matches GATK's FisherStrand: counts summing to more than 2 * min_cell_count are normalised, a total of at
    most min_count gives missing, and p-values are floored at min_p_value. With normalize False, large counts
    get a chi-squared test instead of being normalised.

    Args:
        sb: Either [ref fwd, ref rev, alt fwd, alt rev], or [[ref fwd, ref rev], [alt fwd, alt rev]].
        normalize: Normalise large counts (True), or use a chi-squared test for them (False).
        min_cell_count: The count above which a Fisher exact test isn't used as is.
        min_count: The total count at or below which FS is missing.
        min_p_value: The floor applied to the p-value.

    Returns:
        The FS value.
    """
    if not isinstance(sb, hl.expr.ArrayNumericExpression):
        sb = hl.bind(lambda x: hl.flatten(x), sb)

    sb_sum = hl.bind(lambda x: hl.sum(x), sb)

    # Normalize table if counts get too large
    if normalize:
        fs_expr = hl.bind(
            lambda sb, sb_sum: hl.if_else(
                sb_sum <= 2 * min_cell_count,
                sb,
                sb.map(lambda x: hl.int(x / (sb_sum / min_cell_count))),
            ),
            sb,
            sb_sum,
        )

        # FET
        fs_expr = to_phred(
            hl.max(
                hl.fisher_exact_test(fs_expr[0], fs_expr[1], fs_expr[2], fs_expr[3]).p_value,
                min_p_value,
            )
        )
    else:
        fs_expr = to_phred(
            hl.max(
                hl.contingency_table_test(sb[0], sb[1], sb[2], sb[3], min_cell_count=min_cell_count).p_value,
                min_p_value,
            )
        )

    # Return null if counts <= `min_count`
    return hl.or_missing(
        sb_sum > min_count,
        hl.max(0, fs_expr),  # Needed to avoid -0.0 values
    )


def sor_from_sb(sb: hl.expr.ArrayNumericExpression | hl.expr.ArrayExpression) -> hl.expr.Float64Expression:
    """Computes SOR (symmetric odds ratio) from an SB (strand balance) table, as GATK's StrandOddsRatio does.

    Args:
        sb: Either [ref fwd, ref rev, alt fwd, alt rev], or [[ref fwd, ref rev], [alt fwd, alt rev]].

    Returns:
        The SOR value.
    """
    if not isinstance(sb, hl.expr.ArrayNumericExpression):
        sb = hl.bind(lambda x: hl.flatten(x), sb)

    sb = sb.map(lambda x: hl.float64(x) + 1)

    ref_fw = sb[0]
    ref_rv = sb[1]
    alt_fw = sb[2]
    alt_rv = sb[3]
    symmetrical_ratio = ((ref_fw * alt_rv) / (alt_fw * ref_rv)) + ((alt_fw * ref_rv) / (ref_fw * alt_rv))
    ref_ratio = hl.min(ref_rv, ref_fw) / hl.max(ref_rv, ref_fw)
    alt_ratio = hl.min(alt_fw, alt_rv) / hl.max(alt_fw, alt_rv)
    return hl.log(symmetrical_ratio) + hl.log(ref_ratio) - hl.log(alt_ratio)


def pab_max_expr(
    gt_expr: hl.expr.CallExpression,
    ad_expr: hl.expr.ArrayExpression,
    la_expr: hl.expr.ArrayExpression | None = None,
    n_alleles_expr: hl.expr.Int32Expression | None = None,
) -> hl.expr.ArrayExpression:
    """Aggregates, per alt allele, the maximum p-value of a binomial test of allele balance (PAB) in hets.

    gt_expr and ad_expr may use local or global alleles. For local alleles, pass la_expr and n_alleles_expr
    so they can be converted to global alleles.

    Args:
        gt_expr: The genotype.
        ad_expr: Allele depths.
        la_expr: Local allele indices. If given, gt_expr and ad_expr are converted to global alleles.
        n_alleles_expr: The number of alleles at the site. Required with la_expr.

    Returns:
        The maximum p-value for each alt allele.

    Raises:
        ValueError: If la_expr is given without n_alleles_expr.
    """
    if la_expr is not None:
        if n_alleles_expr is None:
            raise ValueError('Must provide `n_alleles_expr` if `la_expr` is provided!')

        ad_expr = hl.vds.local_to_global(ad_expr, la_expr, n_alleles_expr, fill_value=0, number='R')
        gt_expr = hl.vds.lgt_to_gt(gt_expr, la_expr)

    return hl.agg.array_agg(
        lambda x: hl.agg.filter(
            gt_expr.is_het(),
            hl.agg.max(hl.binom_test(x, hl.sum(ad_expr), 0.5, 'two-sided')),
        ),
        ad_expr[1:],  # Skip ref allele
    )
