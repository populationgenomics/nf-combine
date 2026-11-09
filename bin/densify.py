#!/usr/bin/env python3

"""Hail script, run on the Dataproc cluster by submit_hail.py: one VDS -> a dense MT and two sites-only VCFs.

Ported from cpg-flow-seqr-loader's CreateDenseMtFromVdsWithHail stage (scripts/densify_VDS_to_MT.py):
- the VDS is read on the densify intervals, densified, and filtered to variant sites, then checkpointed
- GATK-style INFO fields are computed and multiallelics are split, each split row keeping only its own alt's
  allele-specific INFO, giving the dense MT
- the MT's rows are exported twice, as a directory of sites-only VCF shards:
  - header per shard: each shard is a valid VCF, so the shards can be processed on their own (VEP, VQSR)
  - separate header: the shards concatenate, after the header, into one VCF (e.g. for VQSR training)

Each export writes a shard-manifest.txt into its directory, listing the shards.

Each step is skipped if its output is complete, so a rerun resumes. The checkpoint and the MT are intermediate
outputs in the temp bucket, kept only for that.
"""

import argparse
import logging

import hail as hl
from cpg_utils import to_path
from nf_combine_hail.intervals import read_bed_file_as_intervals
from nf_combine_hail.sparse_mt import default_compute_info
from nf_combine_hail.vcf import adjust_vcf_incompatible_types

logger = logging.getLogger(__name__)

SHARD_MANIFEST = 'shard-manifest.txt'

# INFO arrays from default_compute_info holding one value per alt allele of the unsplit site
PER_ALT_INFO_FIELDS = [
    'AS_QUALapprox',
    'AS_VarDP',
    'AS_ReadPosRankSum',
    'AS_MQRankSum',
    'AS_MQ',
    'AS_QD',
    'AS_FS',
    'AS_SOR',
    'AS_pab_max',
    'AC',
    'AC_raw',
]
# INFO arrays that describe the whole site, kept as they are on each split row
SITE_ARRAY_INFO_FIELDS = ['SB']


def get_arguments() -> argparse.Namespace:
    """Parses the command line, as built by the DENSIFY process.

    Returns:
        The parsed arguments.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--vds', required=True, help='The VDS to densify')
    parser.add_argument('--intervals', required=True, help='Regions BED to read the VDS on')
    parser.add_argument('--checkpoint', required=True, help='MT checkpoint written after densifying')
    parser.add_argument('--out', required=True, help='The dense MT to write')
    parser.add_argument('--sites_only', required=True, help='Directory of header-per-shard sites-only VCFs')
    parser.add_argument('--separate_header', required=True, help='Directory of separate-header sites-only VCFs')
    return parser.parse_args()


def mt_is_complete(mt_path: str) -> bool:
    """Checks whether a MatrixTable is fully written: Hail writes _SUCCESS at its root last.

    Args:
        mt_path: The MT path.

    Returns:
        True if the MT is complete.
    """
    return to_path(f'{mt_path}/_SUCCESS').exists()


def densify(vds_path: str, intervals_path: str, checkpoint_path: str) -> hl.MatrixTable:
    """Densifies a VDS, keeps its variant sites, and checkpoints the result.

    Args:
        vds_path: The VDS.
        intervals_path: Regions BED, setting the MT's partitions.
        checkpoint_path: Where to checkpoint the dense MT.

    Returns:
        The checkpointed MT.
    """
    vds = hl.vds.read_vds(vds_path, intervals=read_bed_file_as_intervals(intervals_path))
    mt = hl.vds.to_dense_mt(vds)

    # remove monoallelic sites, and sites that are non-ref in no sample
    mt = mt.filter_rows((hl.len(mt.alleles) > 1) & (hl.agg.any(mt.LGT.is_non_ref())))

    # site-level DP, named to avoid a collision with the entry field
    mt = mt.annotate_rows(
        site_dp=hl.agg.sum(mt.DP),
        ANS=hl.agg.count_where(hl.is_defined(mt.LGT)) * 2,
    )
    return mt.checkpoint(checkpoint_path, overwrite=True)


def annotate_and_split(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Adds GATK-style INFO fields to the dense MT, drops the gVCF INFO, and splits multiallelics by allele.

    Args:
        mt: The dense MT, from densify.

    Returns:
        The annotated, split MT.
    """
    info_ht = default_compute_info(mt, site_annotations=True, n_partitions=mt.n_partitions())
    info_ht = info_ht.annotate(info=info_ht.info.annotate(DP=mt.rows()[info_ht.key].site_dp))

    # VQSR can't parse pipe-delimited AS_VarDP (e.g. VarDP=|132.1|140.2): it reads one string value, and fails
    # on the second allele. So no field is pipe-delimited
    info_ht = adjust_vcf_incompatible_types(info_ht, pipe_delimited_annotations=[])

    mt = mt.annotate_rows(info=info_ht[mt.row_key].info)
    mt = mt.drop('gvcf_info')

    logger.info('Splitting multiallelics, in a sparse way')
    return select_split_allele(hl.experimental.sparse_split_multi(mt))


def select_split_allele(mt: hl.MatrixTable) -> hl.MatrixTable:
    """Cuts each split row's allele-specific INFO down to its own alt allele.

    sparse_split_multi copies the unsplit site's INFO onto every split row unchanged, so each row carries the
    values for all of the site's alts. GATK reads an allele-specific value by the alt's position in the record,
    which is always the first in a biallelic row, so without this every split row gets the first alt's values.

    Args:
        mt: The split MT, with a_index, and INFO from default_compute_info after adjust_vcf_incompatible_types.

    Returns:
        The MT, with one value per allele-specific INFO field on each row.

    Raises:
        ValueError: If the INFO has an array field not known to be per-alt or site-level.
    """
    info = mt.info
    arrays = [field for field, dtype in info.dtype.items() if isinstance(dtype, hl.tarray)]
    unknown = sorted(set(arrays) - set(PER_ALT_INFO_FIELDS) - set(SITE_ARRAY_INFO_FIELDS))
    if unknown:
        raise ValueError(f'INFO arrays neither per-alt nor site-level, so not split by allele: {unknown}')

    # a_index is the row's alt in the unsplit site's alleles, so its value is at a_index - 1 among the alts
    selected = {field: info[field][mt.a_index - 1 : mt.a_index] for field in PER_ALT_INFO_FIELDS if field in info}
    if 'AS_SB_TABLE' in info:
        # a string by now, 'ref|alt1|alt2...' (adjust_vcf_incompatible_types): keep the ref and this alt
        tables = info.AS_SB_TABLE.split(r'\|')
        selected['AS_SB_TABLE'] = tables[0] + '|' + tables[mt.a_index]
    return mt.annotate_rows(info=info.annotate(**selected))


def export_sites_only(sites: hl.Table, out_dir: str, parallel: str) -> None:
    """Exports a sites-only VCF as a directory of shards, unless an earlier attempt finished it.

    The shard manifest is written last, so a directory without one is a partial export, which is removed
    before exporting again.

    Args:
        sites: The MT's rows.
        out_dir: The directory to write.
        parallel: Hail's export_vcf parallel mode, 'header_per_shard' or 'separate_header'.
    """
    if to_path(f'{out_dir}/{SHARD_MANIFEST}').exists():
        logger.info(f'Reusing {out_dir}, written by an earlier attempt')
        return
    if to_path(out_dir).exists():
        logger.info(f'Removing the partial export at {out_dir}')
        to_path(out_dir).rmtree()
    logger.info(f'Writing sites-only VCF shards ({parallel}) to {out_dir}')
    hl.export_vcf(sites, out_dir, tabix=True, parallel=parallel)


def main() -> None:
    """Writes the dense MT, then both sites-only VCF exports, reusing whatever an earlier attempt finished."""
    args = get_arguments()
    hl.init(default_reference='GRCh38')

    if mt_is_complete(args.out):
        logger.info(f'Reusing {args.out}, written by an earlier attempt')
    else:
        if mt_is_complete(args.checkpoint):
            logger.info(f'Resuming from the densified checkpoint {args.checkpoint}')
            mt = hl.read_matrix_table(args.checkpoint)
        else:
            logger.info(f'Densifying {args.vds}, checkpointing to {args.checkpoint}')
            mt = densify(args.vds, args.intervals, args.checkpoint)
        logger.info(f'Writing the dense MT to {args.out}')
        annotate_and_split(mt).write(args.out, overwrite=True)

    sites = hl.read_matrix_table(args.out).rows()
    export_sites_only(sites, args.sites_only, parallel='header_per_shard')
    export_sites_only(sites, args.separate_header, parallel='separate_header')


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    main()
