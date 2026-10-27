#!/usr/bin/env python3

"""
Build the combining sample sheet for a set of Metamist cohorts, and find the VDSs already registered in the
project.

Run by the BUILD_SAMPLESHEET process at the start of main.nf, for --cohorts runs only (see DESIGN.md):
- cohorts are looked up within the project; missing, inactive and empty cohorts fail the run
- the cohorts' SGs may belong to other datasets (ADR 0001); each row records the SG's own dataset
- an SG in several cohorts is kept once, and the overlap is logged
- SGs that aren't of the requested sequencing type are skipped with a warning
- each SG needs exactly one distinct active gVCF; several analyses of one file count as one. An SG with none,
  or with several different files, fails the run, and every offender is reported at once (ADR 0002)
- analysis status is ignored: any active gVCF is treated as complete

The sheet is named for the VDS it makes (ADR 0004): <name>.tsv, where the name is the first 38 hex characters
of the SHA-256 of the sorted SG IDs, joined, then 'refresh:' and the sorted refresh IDs if there are any, then
_<SG count>. Without a refresh this is cohort_demo's naming.

Alongside the sheet it writes:
- base_vds.txt: the newest active 'vds' analysis of this sequencing type, or nothing
- registered_vds.txt: every active 'vds' analysis path in the project, one per line

    build_samplesheet.py --project <dataset>[-test] --sequencing-type genome --cohorts COH0001 COH0002 \\
        [--refresh CPG0001]
"""

import csv
import hashlib
import logging
from argparse import ArgumentParser
from pathlib import Path

from metamist.graphql import gql, query

logger = logging.getLogger(__name__)

COLUMNS = ['sg_id', 'gvcf', 'dataset']
SEQUENCING_TYPES = {'genome', 'exome'}

# one round trip: cohort status, each cohort's SGs with their gVCFs, and the project's VDSs
COHORTS_QUERY = """
query BuildCombineSamplesheet($project: String!, $cohorts: [String!]!) {
    project(name: $project) {
        cohorts(id: {in_: $cohorts}) {
            id
            status
            sequencingGroups {
                id
                type
                sample {
                    project {
                        name
                    }
                }
                analyses(type: {eq: "gvcf"}, active: {eq: true}) {
                    id
                    output
                    outputs
                }
            }
        }
        analyses(type: {eq: "vds"}, active: {eq: true}) {
            id
            output
            outputs
            meta
        }
    }
}
"""


class SamplesheetError(Exception):
    """Invalid cohorts, or SGs without exactly one gVCF."""


def analysis_path(analysis: dict) -> str | None:
    """The primary output path: 'output' where Metamist fills it, else the registered outputs block."""
    return analysis.get('output') or (analysis.get('outputs') or {}).get('basename')


def dataset_of(project: str) -> str:
    """The dataset a Metamist project belongs to: its name without any -test suffix."""
    return project.removesuffix('-test')


def check_cohorts(project: str, cohort_ids: list[str], cohorts: list[dict]) -> None:
    """Fails on any cohort that is missing (from this project), inactive, or empty."""
    missing = sorted(set(cohort_ids) - {cohort['id'] for cohort in cohorts})
    inactive = sorted(cohort['id'] for cohort in cohorts if cohort['status'] != 'active')
    empty = sorted(cohort['id'] for cohort in cohorts if not cohort['sequencingGroups'])
    if missing or inactive or empty:
        raise SamplesheetError(
            f'Unusable cohorts in {project} - missing: {missing}, inactive: {inactive}, empty: {empty}',
        )


def build_rows(cohorts: list[dict], sequencing_type: str) -> list[dict]:
    """One row per in-scope SG, each SG once. Reports every SG without exactly one gVCF at once."""
    rows: list[dict] = []
    first_cohort: dict[str, str] = {}
    errors: list[str] = []
    for cohort in cohorts:
        for sg in cohort['sequencingGroups']:
            sg_id = sg['id']
            if sg_id in first_cohort:
                logger.info(f'{sg_id}: in {first_cohort[sg_id]} and {cohort["id"]}, combined once')
                continue
            first_cohort[sg_id] = cohort['id']
            if sg['type'] != sequencing_type:
                logger.warning(f'{sg_id}: skipping {sg["type"]}, this run is {sequencing_type}')
                continue
            gvcfs = sorted({path for analysis in sg['analyses'] if (path := analysis_path(analysis))})
            if len(gvcfs) != 1:
                errors.append(f'{sg_id}: expected one active gVCF, found {len(gvcfs)}: {gvcfs}')
                continue
            rows.append({'sg_id': sg_id, 'gvcf': gvcfs[0], 'dataset': dataset_of(sg['sample']['project']['name'])})
    if errors:
        raise SamplesheetError('\n'.join(errors))
    if not rows:
        raise SamplesheetError(f'No {sequencing_type} sequencing groups to combine')
    return sorted(rows, key=lambda row: row['sg_id'])


def newest_vds(analyses: list[dict], sequencing_type: str) -> str | None:
    """The base VDS for a cohort run: the newest active VDS of this sequencing type, whatever its cohorts."""
    candidates = [
        analysis
        for analysis in analyses
        if analysis_path(analysis) and (analysis.get('meta') or {}).get('sequencing_type') == sequencing_type
    ]
    if not candidates:
        return None
    return analysis_path(max(candidates, key=lambda analysis: analysis['id']))


def registered_vdss(analyses: list[dict]) -> list[str]:
    """Every active VDS path in the project, whatever its sequencing type."""
    return sorted({path for analysis in analyses if (path := analysis_path(analysis))})


def vds_name(sg_ids: list[str], refresh: list[str]) -> str:
    """
    The VDS name for these SGs. A refresh changes a VDS without changing its SGs, so refresh IDs are part of
    the hash; otherwise the refreshed VDS would get its base's name.
    """
    text = ''.join(sorted(sg_ids)) + ('refresh:' + ''.join(sorted(refresh)) if refresh else '')
    return f'{hashlib.sha256(text.encode()).hexdigest()[:38]}_{len(sg_ids)}'


def write_sheet(rows: list[dict], path: str) -> None:
    with Path(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, COLUMNS, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main(project: str, sequencing_type: str, cohort_ids: list[str], refresh: list[str], outdir: str) -> None:
    cohort_ids = sorted(set(cohort_ids))
    result = query(gql(COHORTS_QUERY), variables={'project': project, 'cohorts': cohort_ids})['project']
    check_cohorts(project, cohort_ids, result['cohorts'])
    rows = build_rows(result['cohorts'], sequencing_type)
    base = newest_vds(result['analyses'], sequencing_type)
    registered = registered_vdss(result['analyses'])

    name = vds_name([row['sg_id'] for row in rows], sorted(set(refresh)))
    out = Path(outdir)
    write_sheet(rows, str(out / f'{name}.tsv'))
    (out / 'base_vds.txt').write_text(f'{base}\n' if base else '')
    (out / 'registered_vds.txt').write_text(''.join(f'{path}\n' for path in registered))
    n_datasets = len({row['dataset'] for row in rows})
    logger.info(f'{project} {sequencing_type}: {len(rows)} SGs from {n_datasets} dataset(s) -> {name}, base: {base}')


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--project', required=True, help='Metamist project: the dataset, or <dataset>-test')
    parser.add_argument('--sequencing-type', required=True, choices=sorted(SEQUENCING_TYPES))
    parser.add_argument('--cohorts', required=True, nargs='+', help='Metamist cohort IDs')
    parser.add_argument('--refresh', nargs='*', default=[], help='SG IDs to refresh; they change the VDS name')
    parser.add_argument('--outdir', default='.', help='Where to write the sheet and the VDS lists')
    args = parser.parse_args()
    main(args.project, args.sequencing_type, args.cohorts, args.refresh, args.outdir)
