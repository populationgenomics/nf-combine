"""Metamist-shaped fixtures shared by the sample sheet and stub-run tests. Metamist is never contacted."""

from pathlib import Path

PROJECT = 'ds'


def touch(path: Path) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()
    return str(path)


def gvcf(sg_id: str, analysis_id: int = 1, path: str | None = None) -> dict:
    """An active gVCF analysis, shaped like the sample sheet query's SG analyses."""
    return {'id': analysis_id, 'outputs': {'path': path or f'gs://cpg-ds-main/gvcf/{sg_id}.g.vcf.gz'}}


def sg(sg_id: str, analyses: list[dict] | None = None, project: str = PROJECT, **fields: object) -> dict:
    """A sequencing group, shaped like the sample sheet query's cohort SGs."""
    return {
        'id': sg_id,
        'type': 'genome',
        'sample': {'project': {'name': project}},
        'analyses': [gvcf(sg_id)] if analyses is None else analyses,
    } | fields


def cohort(cohort_id: str, sgs: list[dict], status: str = 'active') -> dict:
    return {'id': cohort_id, 'status': status, 'sequencingGroups': sgs}


def vds(analysis_id: int, path: str, sequencing_type: str = 'genome') -> dict:
    """An active VDS analysis, shaped like the sample sheet query's project analyses."""
    return {'id': analysis_id, 'outputs': {'path': path}, 'meta': {'sequencing_type': sequencing_type}}
