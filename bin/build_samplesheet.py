#!/usr/bin/env python3

"""Build the combining sample sheet for a set of Metamist cohorts, and find the project's registered VDSs.

Run by the BUILD_SAMPLESHEET process at the start of main.nf, for runs built from a list of cohorts:
- cohorts are looked up in Metamist; missing, inactive and empty cohorts fail the run
- Cohort may contain SGs from multiple datasets (ADR 0001); each row records the SG's own dataset for later splitting
- an SG in several cohorts is kept once, and the overlap is logged
- SGs that aren't of the requested sequencing type are skipped with a warning
- each SG needs exactly one distinct active gVCF; several analyses of one file count as one. An SG with none,
  or with several different files, fails the run, and every offender is reported at once (ADR 0002)

The sheet is named for the VDS it makes (ADR 0004): <name>.tsv, where the name is the first 38 hex characters
of the SHA-256 of the sorted SG IDs, joined, then 'refresh:' and the sorted refresh IDs if there are any, then
_<SG count>.

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
from cloudpathlib.anypath import to_anypath
from pathlib import Path

from cpg_utils.existence_checks import exists
from metamist.graphql import gql, query

logger = logging.getLogger(__name__)

COLUMNS = ['sg_id', 'gvcf', 'dataset']
SEQUENCING_TYPES = {'genome', 'exome'}

# one round trip: cohort status, each cohort's SGs with their gVCFs, and the project's VDSs
COHORTS_QUERY = """
query BuildCombineSamplesheet($project: String!, $cohorts: [String!]!) {
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
    project(name: $project) {
        analyses(type: {eq: "combiner"}, active: {eq: true}) {
            id
            outputs
            meta
        }
    }
}
"""


class SamplesheetError(Exception):
    """Raised for invalid cohorts, or for SGs without exactly one gVCF."""


def analysis_path(analysis: dict) -> str | None:
    """Gets the path of a Metamist analysis entry.

    Args:
        analysis: An analysis entry from the GraphQL response.

    Returns:
        The entry's outputs.path, else its output, else None.
    """
    if isinstance(analysis.get('outputs'), str) or analysis.get('outputs') is None:
        return None
    return analysis['outputs'].get('path')


def dataset_of(project: str) -> str:
    """Gets the dataset a Metamist project belongs to.

    Args:
        project: A Metamist project name: the dataset, or <dataset>-test.

    Returns:
        The project name without any -test suffix.
    """
    return project.removesuffix('-test')


def check_cohorts(cohort_ids: list[str], cohorts: list[dict]) -> None:
    """Checks that every requested cohort is usable.

    Args:
        cohort_ids: The requested cohort IDs.
        cohorts: The cohorts the project returned for those IDs.

    Raises:
        SamplesheetError: If any cohort is missing, inactive, or empty.
    """
    missing = sorted(set(cohort_ids) - {cohort['id'] for cohort in cohorts})
    inactive = sorted(cohort['id'] for cohort in cohorts if cohort['status'] != 'active')
    empty = sorted(cohort['id'] for cohort in cohorts if not cohort['sequencingGroups'])
    if missing or inactive or empty:
        raise SamplesheetError(
            f'Unusable cohorts - missing: {missing}, inactive: {inactive}, empty: {empty}',
        )


def build_rows(cohorts: list[dict], sequencing_type: str) -> list[dict]:
    """Builds the sample sheet rows: one per in-scope SG, each SG once.

    An SG in several cohorts is kept once and the overlap logged; SGs of another sequencing type are skipped
    with a warning.

    Args:
        cohorts: Cohorts from the GraphQL response, with their SGs and active gVCFs.
        sequencing_type: The sequencing type being combined, 'genome' or 'exome'.

    Returns:
        Rows with the columns in COLUMNS, sorted by SG ID.

    Raises:
        SamplesheetError: If any SG lacks exactly one distinct active gVCF (every offender is reported at once),
            or if no SGs are left to combine.
    """
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
        raise SamplesheetError(f'No {sequencing_type} sequencing groups to combine.')
    return sorted(rows, key=lambda row: row['sg_id'])


def newest_vds(analyses: list[dict], sequencing_type: str) -> str | None:
    """Finds the base VDS for a cohort run.

    Args:
        analyses: The project's active 'vds' analyses.
        sequencing_type: The sequencing type being combined.

    Returns:
        The path of the newest active VDS of this sequencing type, whatever its cohorts, or None if there is none.
    """
    candidates = [
        analysis
        for analysis in analyses
        if analysis_path(analysis) and (analysis.get('meta') or {}).get('sequencing_type') == sequencing_type
    ]
    if not candidates:
        return None
    return analysis_path(max(candidates, key=lambda analysis: analysis['id']))


def vds_name(sg_ids: list[str], refresh: list[str]) -> str:
    """Names the VDS made from these SGs.

    A refresh changes a VDS without changing its SGs, so refresh IDs are part of the hash; otherwise the
    refreshed VDS would get its base's name.

    Args:
        sg_ids: The SG IDs in the VDS.
        refresh: The SG IDs being refreshed, if any.

    Returns:
        The first 38 hex characters of the hash, then _<SG count>.
    """
    text = ''.join(sorted(sg_ids)) + ('refresh:' + ''.join(sorted(set(refresh))) if refresh else '')
    return f'{hashlib.sha256(text.encode()).hexdigest()[:38]}_{len(sg_ids)}'


def write_sheet(rows: list[dict], path: str | Path) -> None:
    """Writes the sample sheet as a TSV with a header row.

    Args:
        rows: Rows with the columns in COLUMNS.
        path: Where to write the sheet.
    """
    with to_anypath(path).open('w', newline='') as handle:
        writer = csv.DictWriter(handle, COLUMNS, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def validate_file_existence(output_rows: list[dict], base_vds: str | None) -> None:
    """
    Runs real existence checking on the combiner inputs.

    Args:
        output_rows: list of input gVCFs to be combined
        base_vds: nullable, a VDS to build on

    Raises:
        FileNotFoundError: If any expected input file/VDS doesn't exist.
    """

    errors: list[str] = []

    if isinstance(base_vds, str):
        if not exists(base_vds):
            errors.append(f'VDS: {base_vds} does not exist')

    for row in output_rows:
        for ext in ['', '.tbi']:
            filepath = f'{row["gvcf"]}{ext}'
            if not exists(filepath):
                errors.append(f'gVCF: {filepath} does not exist')

    if errors:
        raise FileNotFoundError('\n'.join(errors))


def check_missing_refresh(input_sgids: set[str], refresh_sgids: list[str]) -> None:
    """
    Checks that all SG IDs to be 'refreshed' are present in inputs.
    Args:
        input_sgids: all SG IDs colleted from input Cohorts.
        refresh_sgids: all SG IDs to be refreshed.

    Raises:
        SampleSheetError: If any refresh SG ID is not in the input Cohorts.
    """

    if fresh := (refresh_sgids - input_sgids):
        raise SamplesheetError(f'{len(fresh)} SG IDs were due to be refreshed, but were not in input Cohort(s).')


def main(
    project: str,
    sequencing_type: str,
    cohort_ids: list[str],
    refresh: list[str],
    outdir: str,
    base_vds: str | None = None,
) -> None:
    """Builds the sample sheet, and writes it with the base and registered VDS lists.

    Args:
        project: The Metamist project: the dataset, or <dataset>-test.
        sequencing_type: The sequencing type to combine, 'genome' or 'exome'.
        cohort_ids: The cohorts to combine; duplicates are ignored.
        refresh: SG IDs to refresh; they change the VDS name.
        outdir: Where to write the sheet, base_vds.txt and registered_vds.txt.
        base_vds: Optiona, base VDs to build on. Skips relevant metamist queries.
    """
    cohort_ids = sorted(set(cohort_ids))
    result = query(gql(COHORTS_QUERY), variables={'project': project, 'cohorts': cohort_ids})
    check_cohorts(cohort_ids, result['cohorts'])
    rows = build_rows(result['cohorts'], sequencing_type)

    base = base_vds or newest_vds(result['project']['analyses'], sequencing_type)

    # check for existence of each assay and reference
    validate_file_existence(rows, base)

    # generate a hash from the SG IDs, and optionally IDs to refresh (prevents name collision)
    name = vds_name([row['sg_id'] for row in rows], refresh)

    # check that all SG IDs due to be refreshed are present in inputs
    check_missing_refresh({row['sg_id'] for row in rows}, refresh)

    out = to_anypath(outdir)

    write_sheet(rows, out / f'{name}.tsv')
    (out / 'base_vds.txt').write_text(f'{base}\n' if base else '')
    n_datasets = len({row['dataset'] for row in rows})
    logger.info(f'{project} {sequencing_type}: {len(rows)} SGs from {n_datasets} dataset(s) -> {name}, base: {base}')


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO, format='%(levelname)s %(message)s')
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--project', required=True, help='Metamist project: the dataset, or <dataset>-test')
    parser.add_argument('--sequencing-type', required=True, choices=SEQUENCING_TYPES)
    parser.add_argument('--cohorts', required=True, nargs='+', help='Metamist cohort IDs')
    parser.add_argument('--refresh', nargs='*', default=[], help='SG IDs to refresh; they change the VDS name')
    parser.add_argument('--outdir', default='.', help='Where to write the sheet and the VDS lists')
    parser.add_argument('--base-vds', default=None, help='A base VDS to build on, skips queries.')
    args = parser.parse_args()
    main(args.project, args.sequencing_type, args.cohorts, args.refresh, args.outdir, args.base_vds)
