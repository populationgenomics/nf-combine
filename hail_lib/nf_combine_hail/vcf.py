"""Preparing Hail INFO fields for VCF export.

Copied from cpg-flow-seqr-loader's hail_scripts/vcf.py (originally gnomad_methods), keeping only
adjust_vcf_incompatible_types.
"""

import logging

import hail as hl

logger = logging.getLogger(__name__)

# allele-specific INFO fields GATK expects as pipe-delimited strings, one value per allele
INFO_VCF_AS_PIPE_DELIMITED_FIELDS = [
    'AS_QUALapprox',
    'AS_VarDP',
    'AS_MQ_DP',
    'AS_RAW_MQ',
    'AS_SB_TABLE',
]


def adjust_vcf_incompatible_types(
    ht: hl.Table,
    pipe_delimited_annotations: list[str] = INFO_VCF_AS_PIPE_DELIMITED_FIELDS,  # noqa: B006
) -> hl.Table:
    """Converts a Table's INFO fields to types that can be exported to VCF.

    int64 fields (and arrays of them) are capped and cast to int32. The pipe_delimited_annotations arrays become
    pipe-delimited strings, with an empty leading value for the ref allele as GATK expects. An array-of-arrays SB
    is flattened, and AS_SB_TABLE becomes pipe-delimited comma-joined pairs.

    Args:
        ht: The Table, with an info struct.
        pipe_delimited_annotations: Fields of ht.info to pipe-delimit.

    Returns:
        The Table, keyed by locus and alleles, ready for VCF export.
    """

    def get_pipe_expr(array_expr: hl.expr.ArrayExpression) -> hl.expr.StringExpression:
        return hl.delimit(array_expr.map(lambda x: hl.or_else(hl.str(x), '')), '|')

    # Make sure the HT is keyed by locus, alleles
    ht = ht.key_by('locus', 'alleles')

    info_type_convert_expr = {}
    # Convert int64 fields to int32 (int64 isn't supported by VCF)
    for f, ft in ht.info.dtype.items():
        if ft == hl.dtype('int64'):
            logger.warning(
                'Coercing field info.%s from int64 to int32 for VCF output. Value will be capped at int32 max value.',
                f,
            )
            info_type_convert_expr.update({f: hl.int32(hl.min(2**31 - 1, ht.info[f]))})
        elif ft == hl.dtype('array<int64>'):
            logger.warning(
                'Coercing field info.%s from array<int64> to array<int32> for VCF'
                ' output. Array values will be capped at int32 max value.',
                f,
            )
            info_type_convert_expr.update({f: ht.info[f].map(lambda x: hl.int32(hl.min(2**31 - 1, x)))})

    ht = ht.annotate(info=ht.info.annotate(**info_type_convert_expr))

    info_expr = {}

    # Make sure to pipe-delimit fields that need to.
    # Note: the expr needs to be prefixed by "|" because GATK expect one value for the ref (always empty)
    # Note2: this doesn't produce the correct annotation for AS_SB_TABLE, it
    # is handled below
    for f in pipe_delimited_annotations:
        if f in ht.info and f != 'AS_SB_TABLE':
            info_expr[f] = '|' + get_pipe_expr(ht.info[f])

    # Flatten SB if it is an array of arrays
    if 'SB' in ht.info and not isinstance(ht.info.SB, hl.expr.ArrayNumericExpression):
        info_expr['SB'] = ht.info.SB[0].extend(ht.info.SB[1])

    if 'AS_SB_TABLE' in ht.info:
        info_expr['AS_SB_TABLE'] = get_pipe_expr(ht.info.AS_SB_TABLE.map(lambda x: hl.delimit(x, ',')))

    # Annotate with new expression
    return ht.annotate(info=ht.info.annotate(**info_expr))
