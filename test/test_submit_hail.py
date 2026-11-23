"""bin/submit_hail.py without a cluster: argument splitting, the record, and the completion checks either side of the
job. run_on_cluster is replaced, so Dataproc is never contacted."""

import json
from pathlib import Path

import pytest
import submit_hail
from conftest import touch

CLUSTER_ARGS = [
    '--cluster-name', 'densify-x', '--region', 'r', '--project', 'p', '--staging-bucket', 'gs://s',
    '--temp-bucket', 'gs://t', '--master-type', 'm', '--worker-type', 'w', '--num-workers', '2',
    '--num-secondary-workers', '0', '--max-age', '60', '--max-idle', '60', '--boot-disk-size', '100',
    '--num-local-ssds', '0', '--service-account', 'sa', '--script', 'densify.py',
]  # fmt: skip


def submit(tmp_path: Path, marker: Path, *extra: str) -> dict:
    """Runs main() for one job whose only output is marker; returns the record it wrote."""
    output = tmp_path / 'record.json'
    submit_hail.main(
        [*CLUSTER_ARGS, '--complete', str(marker), '--output', str(output), *extra, '--', '--out', 'x.mt'],
    )
    return json.loads(output.read_text())


def test_split_argv_passes_everything_after_the_first_separator_to_the_script():
    assert submit_hail.split_argv(['--a', '--', '--b', '--', 'c']) == (['--a'], ['--b', '--', 'c'])
    assert submit_hail.split_argv(['--a']) == (['--a'], [])


def test_records_need_key_value_pairs():
    assert submit_hail.parse_records(['vds=gs://b/x.vds', 'base_vds=']) == {'vds': 'gs://b/x.vds', 'base_vds': None}
    with pytest.raises(ValueError, match='key=value'):
        submit_hail.parse_records(['vds'])


def test_complete_outputs_need_no_cluster(tmp_path, monkeypatch):
    def no_cluster(*_args):
        raise AssertionError('a cluster was started')

    monkeypatch.setattr(submit_hail, 'run_on_cluster', no_cluster)
    marker = Path(touch(tmp_path / 'x.mt' / '_SUCCESS'))
    manifest = tmp_path / 'remote-manifest.txt'
    manifest.write_text('part-0\n')
    local = tmp_path / 'sites_only.manifest'

    record = submit(tmp_path, marker, '--fetch', str(manifest), str(local), '--record', 'mt=x.mt')
    assert record == {'mt': 'x.mt', 'cluster': None}
    assert local.read_text() == 'part-0\n'


def test_the_job_runs_with_the_script_arguments(tmp_path, monkeypatch):
    marker = tmp_path / 'x.mt' / '_SUCCESS'
    seen = {}

    def fake_cluster(args, hail_args):
        seen.update(pyfiles=args.pyfiles, hail_args=hail_args)
        touch(marker)
        return {'name': 'densify-x-1234'}

    monkeypatch.setattr(submit_hail, 'run_on_cluster', fake_cluster)
    record = submit(tmp_path, marker, '--pyfiles', 'nf_combine_hail')
    assert record == {'cluster': {'name': 'densify-x-1234'}}
    assert seen == {'pyfiles': ['nf_combine_hail'], 'hail_args': ['--out', 'x.mt']}


def test_a_job_that_leaves_outputs_missing_fails(tmp_path, monkeypatch):
    monkeypatch.setattr(submit_hail, 'run_on_cluster', lambda *_args: {})
    marker = tmp_path / 'x.mt' / '_SUCCESS'
    with pytest.raises(RuntimeError, match='outputs are missing'):
        submit(tmp_path, marker)
    assert not (tmp_path / 'record.json').exists()
