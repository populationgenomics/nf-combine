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
def explicit_vds(tmp_path) -> Path:
    """The --output_vds of a sample-sheet run."""
    return tmp_path / 'elsewhere' / 'one-off.vds'


@pytest.fixture
def run(tmp_path, sheet, explicit_vds):
    """Run main.nf in stub mode, by default for a sample sheet; returns (result, outputDir, trace rows)."""

    def _run(
        *args: str,
        inputs: tuple = ('--samplesheet', str(sheet), '--output_vds', str(explicit_vds)),
        stub: bool = True,
    ):
        out = tmp_path / 'results'
        trace = tmp_path / 'trace.txt'
        trace.unlink(missing_ok=True)
        result = subprocess.run(
            [
                'nextflow', 'run', str(MAIN_NF), *(['-stub'] if stub else []), '-profile', 'test',
                '-with-trace', str(trace),
                '-output-dir', str(out),
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


@pytest.fixture
def cohort_run(tmp_path, sheet, run):
    """
    Run main.nf for --cohorts, with BUILD_SAMPLESHEET's stub copying a sheet (named for its VDS), a base VDS
    and the registered VDSs from a fixtures directory, as if Metamist had returned them.
    """

    def _cohort_run(*args: str, base: str = '', refresh: tuple = ()):
        fixtures = tmp_path / 'fixtures'
        shutil.rmtree(fixtures, ignore_errors=True)
        fixtures.mkdir()
        shutil.copy(sheet, fixtures / f'{vds_name(SGS, list(refresh))}.tsv')
        if refresh:
            args = ('--refresh', ','.join(refresh), *args)
        (fixtures / 'base_vds.txt').write_text(f'{base}\n' if base else '')
        config = tmp_path / 'fixtures.config'
        config.write_text(f"process {{\n    withName: 'BUILD_SAMPLESHEET' {{ ext.fixtures = '{fixtures}' }}\n}}\n")
        return run('-c', str(config), *args, inputs=('--cohorts', 'COH1,COH2'))

    return _cohort_run


def task_names(tasks: list[dict], process: str) -> list[str]:
    return [t['name'] for t in tasks if t['name'].split(' ')[0].endswith(process)]


def task_output(tmp_path: Path, tasks: list[dict], process: str, name: str = '.command.out') -> str:
    workdir = next(t for t in tasks if t['name'].startswith(process))['hash']
    return next((tmp_path / 'work').glob(f'{workdir}*/{name}')).read_text()


def vds_path(out: Path, refresh: list[str] | None = None) -> Path:
    """Where a cohort run writes its VDS: named for its sheet, so for its SGs and refresh IDs."""
    return out / 'vds' / 'genome' / f'{vds_name(SGS, refresh or [])}.vds'


def written_vds(out: Path) -> Path:
    # outputDir is under a symlinked tmp dir on macOS
    return Path(json.loads((out / RUN_DIR / 'output.json').read_text())['vds']).resolve()


def test_samplesheet_run_writes_to_the_explicit_vds(run, explicit_vds):
    result, out, tasks = run()
    assert result.returncode == 0, result.stdout + result.stderr
    assert not task_names(tasks, 'BUILD_SAMPLESHEET')
    assert len(task_names(tasks, 'COMBINE')) == 1
    assert written_vds(out) == explicit_vds.resolve()
    assert json.loads((out / RUN_DIR / 'output.json').read_text())['base_vds'] == ''
    assert (out / RUN_DIR / 'combine_report.json').exists()
    # a sheet passed with --samplesheet isn't republished
    assert not (out / RUN_DIR / 'samplesheet.tsv').exists()


def test_cohort_run_writes_the_vds_its_sheet_is_named_for(cohort_run):
    result, out, _tasks = cohort_run(refresh=('CPG2',))
    assert result.returncode == 0, result.stdout + result.stderr
    assert written_vds(out) == vds_path(out, refresh=['CPG2']).resolve()
    assert (out / RUN_DIR / 'samplesheet.tsv').exists()


def test_refresh_of_an_sg_not_being_combined_fails(run):
    result, _out, tasks = run('--refresh', 'CPG2,CPG9')
    assert result.returncode != 0
    assert "--refresh names SGs that aren't being combined: [CPG9]" in result.stdout + result.stderr
    assert not tasks


def test_missing_gvcfs_and_indexes_are_listed_together(run, tmp_path):
    (tmp_path / 'data' / 'CPG1.g.vcf.gz').unlink()
    (tmp_path / 'data' / 'CPG3.g.vcf.gz.tbi').unlink()
    result, _out, tasks = run()
    console = result.stdout + result.stderr
    assert result.returncode != 0
    assert '2 gVCF file(s) are missing' in console
    assert 'CPG1.g.vcf.gz\n' in console
    assert 'CPG3.g.vcf.gz.tbi' in console
    assert not tasks


@pytest.mark.parametrize(
    ('content', 'message'),
    [
        ('sg_id\tgvcf\nCPG1\tx\n', 'needs the columns [sg_id, gvcf, dataset]'),
        ('sg_id\tgvcf\tdataset\nCPG1\tx\t\n', '1 sample sheet row(s) have empty values'),
        ('sg_id\tgvcf\tdataset\nCPG1\tx\tds\nCPG1\ty\tds\n', 'SGs appear more than once in the sample sheet: [CPG1]'),
        ('sg_id\tgvcf\tdataset\n', 'The sample sheet has no SGs'),
    ],
)
def test_bad_sample_sheets_fail_early(run, sheet, content, message):
    sheet.write_text(content)
    result, _out, tasks = run()
    assert result.returncode != 0
    assert message in result.stdout + result.stderr
    assert not tasks


def test_cohort_run_builds_on_the_discovered_base_and_registers(cohort_run, tmp_path):
    base = tmp_path / 'base.vds'
    base.mkdir()
    result, out, tasks = cohort_run('--register', 'true', base=str(base), refresh=('CPG1',))
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


def test_samplesheet_runs_never_register(run, tmp_path):
    result, _out, tasks = run('--register', 'true')
    assert result.returncode == 0, result.stdout + result.stderr
    assert 'metamist_registration' not in task_output(tmp_path, tasks, 'COMBINE')


def test_force_ignores_the_discovered_base(cohort_run, tmp_path):
    base = tmp_path / 'base.vds'
    base.mkdir()
    result, out, _tasks = cohort_run('--force', 'true', base=str(base))
    assert result.returncode == 0, result.stdout + result.stderr
    assert json.loads((out / RUN_DIR / 'output.json').read_text())['base_vds'] == ''


def test_missing_base_vds_fails(run, tmp_path):
    result, _out, tasks = run('--base_vds', str(tmp_path / 'nowhere.vds'))
    assert result.returncode != 0
    assert "nowhere.vds doesn't exist" in result.stdout + result.stderr
    assert not tasks


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

    result, out, tasks = run(
        '-c', str(config), '--register', 'true', '--refresh', 'CPG3', inputs=('--cohorts', 'COH1,COH2'), stub=False
    )
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
    assert '--intervals gs://cpg-common-main/references/hail_intervals/hg38/' in submit
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
        (('--cohorts', 'COH1'), 'Provide exactly one of --cohorts or --samplesheet'),
        (('--output_vds', 'gs://b/x'), "--output_vds must end in .vds, got 'gs://b/x'"),
        (('--force', 'yes'), "--force must be true or false, got 'yes'"),
        (('--force', 'true', '--base_vds', 'gs://b/x.vds'), '--force makes a fresh build'),
    ],
)
def test_invalid_params_fail_early(run, args, message):
    result, _out, tasks = run(*args)
    assert result.returncode != 0
    assert message in result.stdout + result.stderr
    assert not tasks


@pytest.mark.parametrize(
    'inputs',
    [('--samplesheet', 'sheet.tsv'), ('--cohorts', 'COH1', '--output_vds', 'gs://b/x.vds')],
)
def test_samplesheet_and_output_vds_go_together(run, inputs):
    result, _out, tasks = run(inputs=inputs)
    assert result.returncode != 0
    assert '--samplesheet and --output_vds go together' in result.stdout + result.stderr
    assert not tasks


def test_missing_dataproc_settings_fail_early(run, tmp_path):
    # the test profile sets placeholders, and an empty value can't be passed on the command line
    unset = tmp_path / 'unset.config'
    unset.write_text('params.dataproc_service_account = null\nparams.cpg_container = null\n')
    result, _out, tasks = run('-c', str(unset))
    assert result.returncode != 0
    assert 'Missing required params: --cpg_container, --dataproc_service_account' in result.stdout + result.stderr
    assert not tasks


def test_exome_needs_exome_intervals(run, tmp_path):
    unset = tmp_path / 'unset.config'
    unset.write_text('params.exome_intervals = null\n')
    result, _out, tasks = run('-c', str(unset), '--sequencing_type', 'exome')
    assert result.returncode != 0
    assert 'No intervals for exome: set --exome_intervals' in result.stdout + result.stderr
    assert not tasks
