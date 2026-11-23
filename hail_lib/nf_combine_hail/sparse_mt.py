"""GATK-style INFO fields (site and allele-specific), computed from the gVCF INFO of a dense MT's entries.

Copied from cpg-flow-seqr-loader's hail_scripts/sparse_mt.py (originally gnomad_methods). The logic is
unchanged; the docstrings were rewritten in Google style.
"""

import logging

import hail as hl

from nf_combine_hail.annotations import (
    fs_from_sb,
    get_adj_expr,
    get_lowqual_expr,
    pab_max_expr,
    sor_from_sb,
)

logger = logging.getLogger(__name__)

INFO_AGG_FIELDS = {
    'sum_agg_fields': ['QUALapprox'],
    'int32_sum_agg_fields': ['VarDP'],
    'median_agg_fields': ['ReadPosRankSum', 'MQRankSum'],
    'array_sum_agg_fields': ['SB', 'RAW_MQandDP'],
}

AS_INFO_AGG_FIELDS = {
    'sum_agg_fields': ['AS_QUALapprox', 'AS_RAW_MQ'],
    'int32_sum_agg_fields': ['AS_VarDP'],
    'median_agg_fields': ['AS_RAW_ReadPosRankSum', 'AS_RAW_MQRankSum'],
    'array_sum_agg_fields': ['AS_SB_TABLE'],
}


def _get_info_agg_expr(  # noqa: PLR0915
    mt: hl.MatrixTable,
    sum_agg_fields: list[str] | dict[str, hl.expr.NumericExpression] = INFO_AGG_FIELDS['sum_agg_fields'],
    int32_sum_agg_fields: list[str] | dict[str, hl.expr.NumericExpression] = INFO_AGG_FIELDS['int32_sum_agg_fields'],
    median_agg_fields: list[str] | dict[str, hl.expr.NumericExpression] = INFO_AGG_FIELDS['median_agg_fields'],
    array_sum_agg_fields: list[str] | dict[str, hl.expr.ArrayNumericExpression] = INFO_AGG_FIELDS[
        'array_sum_agg_fields'
    ],
    prefix: str = '',
    treat_fields_as_allele_specific: bool = False,
) -> dict[str, hl.expr.Aggregation]:
    """Creates the aggregators for site-level or allele-specific INFO fields.

    - `SB` in array_sum_agg_fields is aggregated as `AS_SB_TABLE`, following GATK's naming.
    - `RAW_MQandDP` in array_sum_agg_fields, or both `RAW_MQ` and `MQ_DP`, are used to compute `MQ` and then
      dropped, as GATK recommends.
    - Fields given as a list of names must be entry fields of mt, or fields of mt.gvcf_info. An entry field
      wins a name clash.

    Args:
        mt: The input MT.
        sum_agg_fields: Fields to aggregate by sum.
        int32_sum_agg_fields: Fields to aggregate by sum, as int32.
        median_agg_fields: Fields to aggregate by (approximate) median.
        array_sum_agg_fields: Fields to aggregate by element-wise sum over an array.
        prefix: A prefix for the output field names, e.g. 'AS_'.
        treat_fields_as_allele_specific: Treat the fields as having one value per local allele.

    Returns:
        Output field names mapped to their aggregation expressions.

    Raises:
        ValueError: If a named field is in neither the entry schema nor mt.gvcf_info.
    """

    def _agg_list_to_dict(mt: hl.MatrixTable, fields: list[str]) -> dict[str, hl.expr.NumericExpression]:
        out_fields = {}
        if 'gvcf_info' in mt.entry:
            out_fields = {f: mt.gvcf_info[f] for f in fields if f in mt.gvcf_info}

        out_fields.update({f: mt[f] for f in fields if f in mt.entry})

        # Check that all fields were found.
        missing_fields = [f for f in fields if f not in out_fields]
        if missing_fields:
            raise ValueError(
                'Could not find the following field(s)in the MT entry schema (or nested under mt.gvcf_info: {}'.format(
                    ','.join(missing_fields)
                )
            )

        if treat_fields_as_allele_specific:
            # TODO: Change to use hl.vds.local_to_global when fill_value can accept
            #  missing (error in v0.2.119).
            out_fields = {
                f: hl.bind(
                    lambda x: hl.if_else(f == 'AS_SB_TABLE', x, x[1:]),  # noqa: B023
                    hl.range(hl.len(mt.alleles)).map(
                        lambda i: hl.or_missing(mt.LA.contains(i), out_fields[f][mt.LA.index(i)])  # noqa: B023
                    ),
                )
                for f in fields
            }

        return out_fields

    # Map str to expressions where needed.
    if isinstance(sum_agg_fields, list):
        sum_agg_fields = _agg_list_to_dict(mt, sum_agg_fields)

    if isinstance(int32_sum_agg_fields, list):
        int32_sum_agg_fields = _agg_list_to_dict(mt, int32_sum_agg_fields)

    if isinstance(median_agg_fields, list):
        median_agg_fields = _agg_list_to_dict(mt, median_agg_fields)

    if isinstance(array_sum_agg_fields, list):
        array_sum_agg_fields = _agg_list_to_dict(mt, array_sum_agg_fields)

    aggs = [
        (median_agg_fields, lambda x: hl.agg.approx_quantiles(x, 0.5)),
        (sum_agg_fields, hl.agg.sum),
        (int32_sum_agg_fields, lambda x: hl.int32(hl.agg.sum(x))),
        (array_sum_agg_fields, hl.agg.array_sum),
    ]

    # Create aggregators.
    agg_expr = {}
    for agg_fields, agg_func in aggs:
        for k, expr in agg_fields.items():
            if treat_fields_as_allele_specific:
                # If annotation is of the form 'AS_RAW_*_RankSum' it has a histogram
                # representation where keys give the per-variant rank sum value to one
                # decimal place followed by a comma and the corresponding count for
                # that value, so we want to sum the rank sum value (first element).
                # Rename annotation in the form 'AS_RAW_*_RankSum' to 'AS_*_RankSum'.
                if k.startswith('AS_RAW_') and k.endswith('RankSum'):
                    agg_expr[f'{prefix}{k.replace("_RAW", "")}'] = hl.agg.array_agg(
                        lambda x: agg_func(hl.or_missing(hl.is_defined(x), x[0])),  # noqa: B023
                        expr,
                    )
                else:
                    agg_expr[f'{prefix}{k}'] = hl.agg.array_agg(lambda x: agg_func(x), expr)  # noqa: B023
            else:
                agg_expr[f'{prefix}{k}'] = agg_func(expr)

    if treat_fields_as_allele_specific:
        prefix = 'AS_'

    # Handle annotations combinations and casting for specific annotations
    # If RAW_MQandDP is in agg_expr or if both MQ_DP and RAW_MQ are, compute MQ instead
    mq_tuple = None
    if f'{prefix}RAW_MQandDP' in agg_expr:
        logger.info(
            'Computing %sMQ as sqrt(%sRAW_MQandDP[0]/%sRAW_MQandDP[1]). '
            'Note that %sMQ will be set to 0 if %sRAW_MQandDP[1] == 0.',
            *[prefix] * 5,
        )
        mq_tuple = agg_expr.pop(f'{prefix}RAW_MQandDP')
    elif 'AS_RAW_MQ' in agg_expr and treat_fields_as_allele_specific:
        logger.info(
            'Computing AS_MQ as sqrt(AS_RAW_MQ[i]/AD[i+1]). Note that AS_MQ will be set to 0 if AS_RAW_MQ == 0.'
        )
        ad_expr = hl.vds.local_to_global(mt.LAD, mt.LA, hl.len(mt.alleles), fill_value=0, number='R')
        mq_tuple = hl.zip(agg_expr.pop('AS_RAW_MQ'), hl.agg.array_sum(ad_expr[1:]))
    elif f'{prefix}RAW_MQ' in agg_expr and f'{prefix}MQ_DP' in agg_expr:
        logger.info(
            'Computing %sMQ as sqrt(%sRAW_MQ/%sMQ_DP). Note that MQ will be set to 0 if %sRAW_MQ == 0.',
            *[prefix] * 4,
        )
        mq_tuple = (agg_expr.pop(f'{prefix}RAW_MQ'), agg_expr.pop(f'{prefix}MQ_DP'))

    if mq_tuple is not None:
        if treat_fields_as_allele_specific:
            agg_expr[f'{prefix}MQ'] = mq_tuple.map(lambda x: hl.if_else(x[1] > 0, hl.sqrt(x[0] / x[1]), 0))
        else:
            agg_expr[f'{prefix}MQ'] = hl.if_else(mq_tuple[1] > 0, hl.sqrt(mq_tuple[0] / mq_tuple[1]), 0)

    # If both VarDP and QUALapprox are present, also compute QD.
    if f'{prefix}VarDP' in agg_expr and f'{prefix}QUALapprox' in agg_expr:
        logger.info(
            'Computing %sQD as %sQUALapprox/%sVarDP. Note that %sQD will be set to 0 if %sVarDP == 0.',
            *[prefix] * 5,
        )
        var_dp = agg_expr[f'{prefix}VarDP']
        qual_approx = agg_expr[f'{prefix}QUALapprox']
        if treat_fields_as_allele_specific:
            agg_expr[f'{prefix}QD'] = hl.map(
                lambda x: hl.if_else(x[1] > 0, x[0] / x[1], 0),
                hl.zip(qual_approx, var_dp),
            )
        else:
            agg_expr[f'{prefix}QD'] = hl.if_else(var_dp > 0, qual_approx / var_dp, 0)

    # SB needs to be cast to int32 for FS down the line.
    if f'{prefix}SB' in agg_expr:
        agg_expr[f'{prefix}SB'] = agg_expr[f'{prefix}SB'].map(lambda x: hl.int32(x))

    # SB needs to be cast to int32 for FS down the line.
    if 'AS_SB_TABLE' in agg_expr:
        agg_expr['AS_SB_TABLE'] = agg_expr['AS_SB_TABLE'].map(lambda x: x.map(lambda y: hl.int32(y)))

    return agg_expr


def get_as_info_expr(
    mt: hl.MatrixTable,
    sum_agg_fields=INFO_AGG_FIELDS['sum_agg_fields'],
    int32_sum_agg_fields=INFO_AGG_FIELDS['int32_sum_agg_fields'],
    median_agg_fields=INFO_AGG_FIELDS['median_agg_fields'],
    array_sum_agg_fields=INFO_AGG_FIELDS['array_sum_agg_fields'],
    alt_alleles_range_array_field='alt_alleles_range_array',
    treat_fields_as_allele_specific=False,
) -> hl.expr.StructExpression:
    """Aggregates the gVCF INFO fields in a MT's entries into allele-specific VCF INFO fields.

    The field rules are those of _get_info_agg_expr. If treat_fields_as_allele_specific is False, each entry
    field holds one value, which counts towards every global alt allele in the entry's local alleles. If True,
    each entry field holds one value per local allele, mapped to its global allele.

    Args:
        mt: The input MT.
        sum_agg_fields: Fields to aggregate by sum.
        int32_sum_agg_fields: Fields to aggregate by sum, as int32.
        median_agg_fields: Fields to aggregate by (approximate) median.
        array_sum_agg_fields: Fields to aggregate by element-wise sum.
        alt_alleles_range_array_field: A row field holding the alt allele indices,
            e.g. `hl.range(1, hl.len(mt.alleles))`.
        treat_fields_as_allele_specific: Treat the fields as having one value per local allele.

    Returns:
        A struct of the allele-specific INFO fields, plus AS_FS and AS_SOR if there's strand balance data.

    Raises:
        ValueError: If alt_alleles_range_array_field isn't an array<int32> row field.
    """
    if 'DP' in list(sum_agg_fields) + list(int32_sum_agg_fields):
        logger.warning(
            '`DP` was included in allele-specific aggregation, however `DP` is'
            ' typically not aggregated by allele; `VarDP` is.Note that the resulting'
            ' `AS_DP` field will NOT include reference genotypes.'
        )

    agg_expr = _get_info_agg_expr(
        mt=mt,
        sum_agg_fields=sum_agg_fields,
        int32_sum_agg_fields=int32_sum_agg_fields,
        median_agg_fields=median_agg_fields,
        array_sum_agg_fields=array_sum_agg_fields,
        prefix='' if treat_fields_as_allele_specific else 'AS_',
        treat_fields_as_allele_specific=treat_fields_as_allele_specific,
    )

    if alt_alleles_range_array_field not in mt.row or mt[alt_alleles_range_array_field].dtype != hl.dtype(
        'array<int32>'
    ):
        msg = f"'get_as_info_expr' expected a row field '{alt_alleles_range_array_field}' of type array<int32>"
        logger.error(msg)
        raise ValueError(msg)

    if not treat_fields_as_allele_specific:
        # Modify aggregations to aggregate per allele
        agg_expr = {
            f: hl.agg.array_agg(
                lambda ai: hl.agg.filter(mt.LA.contains(ai), expr),  # noqa: B023
                mt[alt_alleles_range_array_field],
            )
            for f, expr in agg_expr.items()
        }

    # Run aggregations
    info = hl.struct(**agg_expr)

    # Add FS and SOR if SB is present.
    if 'AS_SB_TABLE' in info or 'AS_SB' in info:
        drop = []
        # Rename AS_SB to AS_SB_TABLE if present and add SB Ax2 aggregation logic.
        if 'AS_SB' in agg_expr:
            if 'AS_SB_TABLE' in agg_expr:
                logger.warning(
                    'Both `AS_SB` and `AS_SB_TABLE` were specified for aggregation.'
                    ' `AS_SB` will be used for aggregation.'
                )
            as_sb_table = hl.array(
                [
                    info.AS_SB.filter(lambda x: hl.is_defined(x)).fold(lambda i, j: i[:2] + j[:2], [0, 0])  # ref
                ]
            ).extend(
                info.AS_SB.map(lambda x: x[2:])  # each alt
            )
            drop = ['AS_SB']
        else:
            as_sb_table = info.AS_SB_TABLE
        info = info.annotate(
            AS_SB_TABLE=as_sb_table,
            AS_FS=hl.range(1, hl.len(mt.alleles)).map(lambda i: fs_from_sb(as_sb_table[0].extend(as_sb_table[i]))),
            AS_SOR=hl.range(1, hl.len(mt.alleles)).map(lambda i: sor_from_sb(as_sb_table[0].extend(as_sb_table[i]))),
        ).drop(*drop)

    return info


def get_site_info_expr(
    mt: hl.MatrixTable,
    sum_agg_fields: list[str] | dict[str, hl.expr.NumericExpression] = INFO_AGG_FIELDS['sum_agg_fields'],
    int32_sum_agg_fields: list[str] | dict[str, hl.expr.NumericExpression] = INFO_AGG_FIELDS['int32_sum_agg_fields'],
    median_agg_fields: list[str] | dict[str, hl.expr.NumericExpression] = INFO_AGG_FIELDS['median_agg_fields'],
    array_sum_agg_fields: list[str] | dict[str, hl.expr.ArrayNumericExpression] = INFO_AGG_FIELDS[
        'array_sum_agg_fields'
    ],
) -> hl.expr.StructExpression:
    """Aggregates the gVCF INFO fields in a MT's entries into site-level VCF INFO fields.

    The field rules are those of _get_info_agg_expr. Fields are aggregated over non-ref genotypes only, except
    DP, which needs a dense MT and is aggregated over all genotypes.

    Args:
        mt: The input MT.
        sum_agg_fields: Fields to aggregate by sum.
        int32_sum_agg_fields: Fields to aggregate by sum, as int32.
        median_agg_fields: Fields to aggregate by (approximate) median.
        array_sum_agg_fields: Fields to aggregate by element-wise sum.

    Returns:
        A struct of the site-level INFO fields, plus FS and SOR if SB is aggregated.
    """
    if 'DP' in list(sum_agg_fields) + list(int32_sum_agg_fields):
        logger.warning(
            '`DP` was included in site-level aggregation. This requires a densifying'
            ' prior to running get_site_info_expr'
        )

    agg_expr = _get_info_agg_expr(
        mt=mt,
        sum_agg_fields=sum_agg_fields,
        int32_sum_agg_fields=int32_sum_agg_fields,
        median_agg_fields=median_agg_fields,
        array_sum_agg_fields=array_sum_agg_fields,
    )

    # Add FS and SOR if SB is present
    # This is done outside _get_info_agg_expr as the behavior is different
    # in site vs allele-specific versions
    if 'SB' in agg_expr:
        agg_expr['FS'] = fs_from_sb(agg_expr['SB'])
        agg_expr['SOR'] = sor_from_sb(agg_expr['SB'])

    # Run aggregator on non-ref genotypes
    info = hl.agg.filter(
        mt.LGT.is_non_ref(),
        hl.struct(**{k: v for k, v in agg_expr.items() if k != 'DP'}),
    )

    # Add DP, computed over both ref and non-ref genotypes, if present
    if 'DP' in agg_expr:
        info = info.annotate(DP=agg_expr['DP'])

    return info


def default_compute_info(
    mt: hl.MatrixTable,
    site_annotations: bool = False,
    as_annotations: bool = False,
    # Set to True by default to prevent a breaking change.
    quasi_as_annotations: bool = True,
    n_partitions: int | None = 5000,
    lowqual_indel_phred_het_prior: int = 40,
    ac_filter_groups: dict[str, hl.Expression] | None = None,
) -> hl.Table:
    """Computes a Table of the usual GATK allele-specific INFO fields, with ACs and lowqual flags.

    Multiallelic sites are not split.

    Args:
        mt: The input MT, filtered to sites that are non-ref in some sample.
        site_annotations: Also compute site-level INFO fields.
        as_annotations: Compute allele-specific fields from the allele-specific annotations in gvcf_info.
        quasi_as_annotations: Compute allele-specific fields from the non-allele-specific annotations in
            gvcf_info, aggregated per allele. This approximates allele-specific annotations where the genotypes
            lack them.
        n_partitions: If set, the output is naive-coalesced to this many partitions.
        lowqual_indel_phred_het_prior: Phred-scaled prior for a het genotype at a low-quality indel site. 40
            (1/10k bases) matches the Broad Data Sciences Platform's filtering for VQSR.
        ac_filter_groups: Sample filter expressions, each giving an extra grouping of ACs.

    Returns:
        A Table of the INFO fields, keyed by locus and alleles.

    Raises:
        ValueError: If site_annotations, as_annotations and quasi_as_annotations are all False.
    """
    if not site_annotations and not as_annotations and not quasi_as_annotations:
        raise ValueError(
            'At least one of `site_annotations`, `as_annotations`, or `quasi_as_annotations` must be True!'
        )

    # Add a temporary annotation for allele count groupings.
    ac_filter_groups = {'': True, **(ac_filter_groups or {})}
    mt = mt.annotate_cols(_ac_filter_groups=ac_filter_groups)

    # Move gvcf info entries out from nested struct.
    mt = mt.transmute_entries(**mt.gvcf_info)

    # Adding alt_alleles_range_array as a required annotation for
    # get_as_info_expr to reduce memory usage.
    mt = mt.annotate_rows(alt_alleles_range_array=hl.range(1, hl.len(mt.alleles)))

    info_expr = None
    quasi_info_expr = None

    # Compute quasi-AS info expr.
    if quasi_as_annotations:
        info_expr = get_as_info_expr(mt)

    # Compute AS info expr using gvcf_info allele specific annotations.
    if as_annotations:
        if info_expr is not None:
            quasi_info_expr = info_expr
        info_expr = get_as_info_expr(
            mt,
            **AS_INFO_AGG_FIELDS,
            treat_fields_as_allele_specific=True,
        )

    if info_expr is not None:
        # Add allele specific pab_max
        info_expr = info_expr.annotate(AS_pab_max=pab_max_expr(mt.LGT, mt.LAD, mt.LA, hl.len(mt.alleles)))

    if site_annotations:
        site_expr = get_site_info_expr(mt)
        info_expr = site_expr if info_expr is None else info_expr.annotate(**site_expr)

    # Add 'AC' and 'AC_raw' for each allele count filter group requested.
    # First compute ACs for each non-ref allele, grouped by adj.
    grp_ac_expr = {
        f: hl.agg.array_agg(
            lambda ai: hl.agg.filter(
                mt.LA.contains(ai) & mt._ac_filter_groups[f],  # noqa: B023,SLF001
                hl.agg.group_by(
                    get_adj_expr(mt.LGT, mt.GQ, mt.DP, mt.LAD),
                    hl.agg.sum(mt.LGT.one_hot_alleles(mt.LA.map(lambda x: hl.str(x)))[mt.LA.index(ai)]),
                ),
            ),
            mt.alt_alleles_range_array,
        )
        for f in ac_filter_groups
    }

    # Then, for each non-ref allele, compute
    # 'AC' as the adj group
    # 'AC_raw' as the sum of adj and non-adj groups
    info_expr = info_expr.annotate(
        **{
            f'AC{"_" + f if f else f}_raw': grp.map(lambda i: hl.int32(i.get(True, 0) + i.get(False, 0)))
            for f, grp in grp_ac_expr.items()
        },
        **{f'AC{"_" + f if f else f}': grp.map(lambda i: hl.int32(i.get(True, 0))) for f, grp in grp_ac_expr.items()},
    )

    ann_expr = {'info': info_expr}
    if quasi_info_expr is not None:
        ann_expr['quasi_info'] = quasi_info_expr

    info_ht = mt.select_rows(**ann_expr).rows()

    # Add AS lowqual flag
    info_ht = info_ht.annotate(
        AS_lowqual=get_lowqual_expr(
            info_ht.alleles,
            info_ht.info.AS_QUALapprox,
            indel_phred_het_prior=lowqual_indel_phred_het_prior,
        )
    )

    if site_annotations:
        # Add lowqual flag
        info_ht = info_ht.annotate(
            lowqual=get_lowqual_expr(
                info_ht.alleles,
                info_ht.info.QUALapprox,
                indel_phred_het_prior=lowqual_indel_phred_het_prior,
            )
        )

    if n_partitions is not None:
        info_ht = info_ht.naive_coalesce(n_partitions)

    return info_ht
