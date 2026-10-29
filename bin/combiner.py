#!/usr/bin/env python3

"""Hail script, run on the Dataproc cluster by submit_combiner.py: gVCFs, and optionally a base VDS, -> one VDS.

With a base VDS (see DESIGN.md):
- SGs in the base that aren't in the sample sheet are removed
- refreshed SGs that are in the base are removed, then added back from their current gVCF
- sample sheet SGs that the base lacks are added from their gVCFs

Writes a JSON report of the base VDS used and the SGs added, removed and refreshed. The combiner plan and
temp files are keyed by the VDS name, so a rerun with the same inputs resumes rather than starting again.
"""

import argparse
import gzip
import json
import logging
import os
from csv import DictReader

import hail as hl
from cpg_utils import to_path
from hail.vds import new_combiner

logger = logging.getLogger(__name__)


def get_arguments() -> argparse.Namespace:
    """Parses the command line, as built by submit_combiner.py.

    Returns:
        The parsed arguments.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samplesheet', required=True)
    parser.add_argument('--out', required=True)
    parser.add_argument('--temp_dir', required=True)
    parser.add_argument('--report', required=True)
    parser.add_argument('--intervals', required=True)
    parser.add_argument('--base_vds', default=None)
    parser.add_argument('--refresh', nargs='*', default=[])
    parser.add_argument('--force', action='store_true')
    parser.add_argument('--branch', type=int, required=True)
    parser.add_argument('--gvcf_batch', type=int, required=True)
    parser.add_argument('--target_records', type=int, required=True)
    return parser.parse_args()


def parse_samplesheet(ss_path: str) -> dict[str, str]:
    """Reads the sample sheet.

    Args:
        ss_path: Path to the sample sheet TSV, with sg_id and gvcf columns.

    Returns:
        A mapping of SG ID to gVCF path.
    """
    with open(ss_path) as handle:
        return {row['sg_id']: row['gvcf'] for row in DictReader(handle, delimiter='\t')}


def read_bed_file_as_intervals(bed_path: str) -> list[hl.Interval]:
    """Reads a BED file as GRCh38 intervals.

    Args:
        bed_path: Path to a plain-text or gzipped BED file. Blank lines and '#' lines are skipped.

    Returns:
        One closed interval per BED record, converted from BED's 0-based, half-open coordinates.
    """
    path = to_path(bed_path)
    intervals: list[hl.Interval] = []
    with gzip.open(path.open('rb'), 'rt') if str(path).endswith('.gz') else path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            chrom, start, end = stripped.split()[:3]
            start_locus = hl.Locus(chrom, int(start) + 1, reference_genome='GRCh38')
            end_locus = hl.Locus(chrom, int(end), reference_genome='GRCh38')
            intervals.append(hl.Interval(start_locus, end_locus, includes_start=True, includes_end=True))
    return intervals


def vds_is_complete(vds_path: str) -> bool:
    """Checks whether a VDS is fully written.

    Hail marks each of a VDS's two matrix tables with _SUCCESS once it's fully written.

    Args:
        vds_path: The VDS path.

    Returns:
        True if both matrix tables are marked complete.
    """
    return all(to_path(f'{vds_path}/{part}/_SUCCESS').exists() for part in ['reference_data', 'variant_data'])


def remove_sgs(vds_path: str, sg_ids: set[str]) -> hl.vds.VariantDataset:
    """Reads the VDS without these SGs. Needed for withdrawn and refreshed SGs.

    Args:
        vds_path: The base VDS.
        sg_ids: The SGs to remove.

    Returns:
        The VDS without these SGs, and without any alleles only they carried.
    """
    logger.info(f'Removing {len(sg_ids)} SG(s) from {vds_path}: {sorted(sg_ids)}')
    vds = hl.vds.read_vds(vds_path)
    return hl.vds.filter_samples(vds, sorted(sg_ids), keep=False, remove_dead_alleles=True)


def main() -> None:
    """Combines the sample sheet's gVCFs, and any base VDS, into one VDS, then writes the report.

    Raises:
        RuntimeError: If the output VDS already exists; it can only be a partial write.
    """
    args = get_arguments()
    hl.init(default_reference='GRCh38')

    # submit_combiner.py skips the cluster for a complete VDS, so anything here is a partial write
    if to_path(args.out).exists():
        raise RuntimeError(f'{args.out} exists but is incomplete. Inspect and delete it, then rerun')

    name = to_path(args.out).stem
    combiner_temp = os.path.join(args.temp_dir, 'combiner_temp', name)
    plan = os.path.join(args.temp_dir, 'combiner_plan', f'{name}.json')

    sg_to_gvcf = parse_samplesheet(args.samplesheet)
    requested = set(sg_to_gvcf)
    refresh = set(args.refresh)

    base_sgs: set[str] = set()
    if args.base_vds:
        base_sgs = set(hl.vds.read_vds(args.base_vds).variant_data.s.collect())
    refreshed = refresh & base_sgs
    if refresh - base_sgs:
        logger.warning(f'Refreshed SGs not in the base VDS are simply added: {sorted(refresh - base_sgs)}')
    removed = (base_sgs - requested) | refreshed
    added = sorted(requested - (base_sgs - removed))

    force = args.force
    vds_paths = [args.base_vds] if args.base_vds else []
    if removed:
        # the base minus the removed SGs, kept so a retry can skip the removal and resume the combiner plan
        removal_vds = os.path.join(args.temp_dir, 'combiner_removal', f'{name}.vds')
        if not added:
            # nothing to combine, so the removal is the whole job
            remove_sgs(args.base_vds, removed).write(args.out)
        elif vds_is_complete(removal_vds):
            logger.info(f'Reusing {removal_vds}, written by an earlier attempt')
            vds_paths = [removal_vds]
        else:
            remove_sgs(args.base_vds, removed).write(removal_vds, overwrite=True)
            vds_paths = [removal_vds]
            # the combiner's base has changed, so any saved plan is stale
            force = True

    if added:
        gvcfs = [sg_to_gvcf[sg_id] for sg_id in added]
        # new_combiner resumes from save_path itself, reapplying the tuning args over the saved plan
        combiner = new_combiner(
            output_path=args.out,
            save_path=plan,
            gvcf_paths=gvcfs,
            gvcf_sample_names=added,
            gvcf_external_header=gvcfs[0],
            vds_paths=vds_paths,
            reference_genome='GRCh38',
            temp_path=combiner_temp,
            intervals=read_bed_file_as_intervals(args.intervals),
            force=force,
            branch_factor=args.branch,
            target_records=args.target_records,
            gvcf_batch_size=args.gvcf_batch,
        )
        combiner.run()
    elif not removed:
        # an explicit base holding exactly the requested SGs
        # shouldn't be possible, consider making this a failing case (equivalent SGs should give the same path)
        hl.vds.read_vds(args.base_vds).write(args.out)

    report = {
        'vds': args.out,
        'base_vds': args.base_vds,
        'added': added,
        'removed': sorted(removed - refreshed),
        'refreshed': sorted(refreshed),
        'n_sgs': len(requested),
    }
    with to_path(args.report).open('w') as handle:
        json.dump(report, handle, indent=2)
    logger.info(f'{args.out}: {len(added)} added, {len(report["removed"])} removed, {len(refreshed)} refreshed')


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    main()
