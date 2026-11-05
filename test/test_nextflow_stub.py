"""Offline `nextflow run` checks of main.nf's wiring: the pre-flight checks, where the VDS is written, base VDS
choice, the no-op for an already-combined SG set, registration and the run records. Placeholder files stand in
for gVCFs. Most runs are stubs; one real run swaps in fake executables to check the commands as rendered. The
tests are skipped when nextflow isn't on PATH."""

import csv
import json
import shutil
import subprocess
from pathlib import Path

import pytest
from build_samplesheet import vds_name, write_sheet
from conftest import PROJECT, touch

pytestmark = pytest.mark.skipif(shutil.which('nextflow') is None, reason='nextflow not installed')

MAIN_NF = Path(__file__).parents[1] / 'main.nf'
SGS = ['CPG1', 'CPG2', 'CPG3']
RUN_DIR = Path('combiner_samplesheets/local-test')


def read_tsv(path: Path) -> list[dict]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


@pytest.fixture
def sheet(tmp_path) -> Path:
    """A sample sheet of local placeholder gVCFs, each with its index."""
    data = tmp_path / 'data'
    rows = []
    for sg_id in SGS:
        touch(data / f'{sg_id}.g.vcf.gz.tbi')
        rows.append({'sg_id': sg_id, 'gvcf': touch(data / f'{sg_id}.g.vcf.gz'), 'dataset': PROJECT})
    path = tmp_path / 'sheet.tsv'
    write_sheet(rows, str(path))
    return path


@pytest.fixture
def run(tmp_path, sheet):
    """
    Run main.nf for --cohorts, by default in stub mode; returns (result, outputDir, trace rows).
    BUILD_SAMPLESHEET's stub copies a sheet (named for its VDS), and a base VDS, from a fixtures directory, as
    if Metamist had returned them.
    """

    def _run(
        *args: str,
        base: str = '',
        refresh: tuple = (),
        inputs: tuple = ('--cohorts', 'COH1,COH2'),
        stub: bool = True,
    ):
        fixtures = tmp_path / 'fixtures'
        shutil.rmtree(fixtures, ignore_errors=True)
        fixtures.mkdir()
        shutil.copy(sheet, fixtures / f'{vds_name(SGS, list(refresh))}.tsv')
        if refresh:
            args = ('--refresh', ','.join(refresh), *args)
        (fixtures / 'base_vds.txt').write_text(f'{base}\n' if base else '')
        config = tmp_path / 'fixtures.config'
        config.write_text(f"process {{\n    withName: 'BUILD_SAMPLESHEET' {{ ext.fixtures = '{fixtures}' }}\n}}\n")

        out = tmp_path / 'results'
        trace = tmp_path / 'trace.txt'
        trace.unlink(missing_ok=True)
        result = subprocess.run(
            [
                'nextflow', 'run', str(MAIN_NF), *(['-stub'] if stub else []), '-profile', 'test',
                '-with-trace', str(trace),
                '-output-dir', str(out),
                '-c', str(config),
                *inputs,
                '--dataset', PROJECT,
                '--sequencing_type', 'genome',
                *args,
            ],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            check=False,
        )  # fmt: skip
        tasks = read_tsv(trace) if trace.exists() else []
        return result, out, tasks

    return _run


def task_names(tasks: list[dict], process: str) -> list[str]:
    return [t['name'] for t in tasks if t['name'].split(' ')[0].endswith(process)]


def task_output(tmp_path: Path, tasks: list[dict], process: str, name: str = '.command.out') -> str:
    workdir = next(t for t in tasks if t['name'].startswith(process))['hash']
    return next((tmp_path / 'work').glob(f'{workdir}*/{name}')).read_text()


def vds_path(out: Path, refresh: list[str] | None = None) -> Path:
    """Where a run writes its VDS: named for its sheet, so for its SGs and refresh IDs."""
    return out / 'vds' / 'genome' / f'{vds_name(SGS, refresh or [])}.vds'


def written_vds(out: Path) -> Path:
    # outputDir is under a symlinked tmp dir on macOS
    return Path(json.loads((out / RUN_DIR / 'output.json').read_text())['vds']).resolve()


def test_run_writes_the_vds_its_sheet_is_named_for(run):
    result, out, tasks = run(refresh=('CPG2',))
    assert result.returncode == 0, result.stdout + result.stderr
    assert len(task_names(tasks, 'COMBINE')) == 1
    assert written_vds(out) == vds_path(out, refresh=['CPG2']).resolve()
    assert json.loads((out / RUN_DIR / 'output.json').read_text())['base_vds'] == ''
    assert (out / RUN_DIR / 'combine_report.json').exists()
    assert (out / RUN_DIR / 'samplesheet.tsv').exists()


def test_run_builds_on_the_discovered_base_and_registers(run, tmp_path):
    base = tmp_path / 'base.vds'
    base.mkdir()
    result, out, tasks = run('--register', 'true', base=str(base), refresh=('CPG1',))
    assert result.returncode == 0, result.stdout + result.stderr
    assert len(task_names(tasks, 'BUILD_SAMPLESHEET')) == 1

    cluster = json.loads((out / RUN_DIR / 'output.json').read_text())
    assert cluster['base_vds'] == str(base)
    register = task_output(tmp_path, tasks, 'COMBINE')
    assert '--dry --project ds ' in register
    assert '--type vds --cohorts COH1 COH2 ' in register
    assert 'sequencing_type=genome ar-guid=local-test' in register
    assert f'base_vds={base}' in register
    assert 'refresh=CPG1' in register


def test_register_false_skips_registration(run, tmp_path):
    result, _out, tasks = run('--register', 'false')
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'metamist_registration' not in task_output(tmp_path, tasks, 'COMBINE')


def test_force_ignores_the_discovered_base(run, tmp_path):
    base = tmp_path / 'base.vds'
    base.mkdir()
    result, out, _tasks = run('--force', 'true', base=str(base))
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads((out / RUN_DIR / 'output.json').read_text())['base_vds'] == ''


def test_real_commands_as_rendered(run, sheet, tmp_path):
    """
    A real (not stub) cohort run, with build_samplesheet.py, submit_combiner.py and python3 replaced by fakes
    on PATH: the sheet is copied in, and the submitter and registration only record their arguments.
    """
    base = tmp_path / 'base.vds'
    base.mkdir()
    fake = tmp_path / 'fake'
    fakes = {
        'build_samplesheet.py': (
            f'cp {sheet} {vds_name(SGS, ["CPG3"])}.tsv; echo {base} > base_vds.txt; touch registered_vds.txt'
        ),
        'submit_combiner.py': 'echo {} > output.json; echo {} > combine_report.json',
        'python3': '',
    }
    for name, body in fakes.items():
        touch(fake / name)
        (fake / name).write_text(f'#!/bin/bash\necho "$@" > {name}.args\n{body}\n')
        (fake / name).chmod(0o755)
    config = tmp_path / 'fake.config'
    config.write_text(f"process.beforeScript = 'export PATH={fake}:$PATH'\n")

    result, out, tasks = run('-c', str(config), '--register', 'true', '--refresh', 'CPG3', stub=False)
    assert result.returncode == 0, result.stdout + result.stderr

    build = task_output(tmp_path, tasks, 'BUILD_SAMPLESHEET', 'build_samplesheet.py.args')
    assert build.split() == [
        '--project', 'ds', '--sequencing-type', 'genome', '--cohorts', 'COH1', 'COH2', '--refresh', 'CPG3',
    ]  # fmt: skip

    submit = task_output(tmp_path, tasks, 'COMBINE', 'submit_combiner.py.args')
    assert '--cluster-name combine-local-test ' in submit
    assert '--script /' in submit  # the Hail script was found on PATH
    words = submit.split()
    assert Path(words[words.index('--vds') + 1]).resolve() == vds_path(out, refresh=['CPG3']).resolve()
    assert f'--base-vds {base} ' in submit
    assert '--refresh CPG3 ' in submit
    assert '--temp-bucket gs://cpg-ds-main-tmp ' in submit
    assert '--intervals gs://cpg-common-main/rd_references/genome_combiner_intervals/' in submit
    assert '--preemptible-secondary-workers ' in submit
    assert '--force' not in submit
    assert submit.strip().endswith('--labels ar-guid=local-test compute-category=dataproc dataset=ds')

    register = task_output(tmp_path, tasks, 'COMBINE', 'python3.args')
    assert register.startswith('-m cpg_utils.metamist_registration --project ds --output ')
    assert '--type vds --cohorts COH1 COH2 ' in register
    assert register.strip().endswith(f'--meta sequencing_type=genome ar-guid=local-test base_vds={base} refresh=CPG3')


@pytest.mark.parametrize(
    ('args', 'message'),
    [
        (('--sequencing_type', 'wgs'), "--sequencing_type must be genome or exome, got 'wgs'"),
        (('--accessLevel', 'prod'), "--accessLevel must be full, standard or test, got 'prod'"),
        (('--force', 'yes'), "--force must be true or false, got 'yes'"),
        (('--force', 'true', '--base_vds', 'gs://b/x.vds'), '--force makes a fresh build'),
    ],
)
def test_invalid_params_fail_early(run, args, message):
    result, _out, tasks = run(*args)
    assert result.returncode != 0
    assert message in result.stdout + result.stderr
    assert not tasks


def test_cohorts_are_required(run):
    result, _out, tasks = run(inputs=())
    assert result.returncode != 0
    assert 'Missing required params: --cohorts' in result.stdout + result.stderr
    assert not tasks


def test_missing_dataproc_settings_fail_early(run, tmp_path):
    # the test profile sets placeholders, and an empty value can't be passed on the command line
    unset = tmp_path / 'unset.config'
    unset.write_text('params.cpg_container = null\n')
    result, _out, tasks = run('-c', str(unset))
    assert result.returncode != 0
    assert 'Missing required params: --cpg_container' in result.stdout + result.stderr
    assert not tasks


def test_exome_needs_exome_combiner_intervals(run, tmp_path):
    unset = tmp_path / 'unset.config'
    unset.write_text('params.exome_combiner_intervals = null\n')
    result, _out, tasks = run('-c', str(unset), '--sequencing_type', 'exome')
    assert result.returncode != 0
    assert 'No combiner intervals set for exome: set --exome_combiner_intervals' in result.stdout + result.stderr
    assert not tasks
