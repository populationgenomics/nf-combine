"""bin/build_samplesheet.py: Metamist responses -> sample sheet rows and the project's VDSs. The query is validated
against the schema bundled with the metamist client; responses are hand-built dicts in the shape it returns."""

import csv
import logging

import build_samplesheet
import pytest
from build_samplesheet import (
    COHORTS_QUERY,
    SamplesheetError,
    analysis_path,
    build_rows,
    check_cohorts,
    newest_vds,
    validate_file_existence,
    vds_name,
)
from conftest import PROJECT, cohort, gvcf, sg, vds
from metamist.graphql import gql, validate


@pytest.fixture(autouse=True)
def missing(monkeypatch) -> set[str]:
    """Paths the patched existence check reports as missing; every other path exists. Storage is never contacted."""
    absent: set[str] = set()
    monkeypatch.setattr(build_samplesheet, 'exists', lambda path: path not in absent)
    return absent


def test_query_matches_metamist_schema():
    assert validate(gql(COHORTS_QUERY), use_local_schema=True)


def test_one_row_per_sg_sorted_with_its_own_dataset():
    cohorts = [cohort('COH1', [sg('CPG3'), sg('CPG1', project='other-test')]), cohort('COH2', [sg('CPG2')])]
    assert build_rows(cohorts, 'genome') == [
        {'sg_id': 'CPG1', 'gvcf': 'gs://cpg-ds-main/gvcf/CPG1.g.vcf.gz', 'dataset': 'other'},
        {'sg_id': 'CPG2', 'gvcf': 'gs://cpg-ds-main/gvcf/CPG2.g.vcf.gz', 'dataset': PROJECT},
        {'sg_id': 'CPG3', 'gvcf': 'gs://cpg-ds-main/gvcf/CPG3.g.vcf.gz', 'dataset': PROJECT},
    ]


def test_sg_in_two_cohorts_is_combined_once_and_logged(caplog):
    caplog.set_level(logging.INFO)
    rows = build_rows([cohort('COH1', [sg('CPG1')]), cohort('COH2', [sg('CPG1'), sg('CPG2')])], 'genome')
    assert [row['sg_id'] for row in rows] == ['CPG1', 'CPG2']
    assert 'CPG1: in COH1 and COH2, combined once' in caplog.text


def test_other_sequencing_types_are_skipped_with_a_warning(caplog):
    rows = build_rows([cohort('COH1', [sg('CPG1'), sg('CPG2', type='exome')])], 'genome')
    assert [row['sg_id'] for row in rows] == ['CPG1']
    assert 'CPG2: skipping exome' in caplog.text


def test_several_registrations_of_one_gvcf_count_as_one():
    analyses = [gvcf('CPG1', 1), gvcf('CPG1', 2)]
    assert build_rows([cohort('COH1', [sg('CPG1', analyses)])], 'genome')[0]['gvcf'].endswith('CPG1.g.vcf.gz')


def test_every_sg_without_exactly_one_gvcf_is_reported_at_once():
    sgs = [
        sg('CPG1', []),
        sg('CPG2', [gvcf('CPG2', 1), gvcf('CPG2', 2, 'gs://b/CPG2.recalled.g.vcf.gz')]),
        sg('CPG3'),
    ]
    with pytest.raises(SamplesheetError) as error:
        build_rows([cohort('COH1', sgs)], 'genome')
    assert 'CPG1: expected one active gVCF, found 0' in str(error.value)
    assert 'CPG2: expected one active gVCF, found 2' in str(error.value)
    assert 'CPG3' not in str(error.value)


def test_no_sgs_of_the_sequencing_type_fails():
    with pytest.raises(SamplesheetError, match='No exome sequencing groups to combine'):
        build_rows([cohort('COH1', [sg('CPG1')])], 'exome')


def test_unusable_cohorts_are_reported_together():
    cohorts = [cohort('COH2', [sg('CPG1')], status='inactive'), cohort('COH3', [])]
    with pytest.raises(SamplesheetError) as error:
        check_cohorts(['COH1', 'COH2', 'COH3'], cohorts)
    assert "missing: ['COH1'], inactive: ['COH2'], empty: ['COH3']" in str(error.value)


def test_analysis_path_falls_back_to_the_outputs_block():
    assert analysis_path({'output': None, 'outputs': {'path': 'gs://b/x.vds'}}) == 'gs://b/x.vds'
    assert analysis_path({'output': None, 'outputs': None}) is None


def test_base_is_the_newest_vds_of_the_sequencing_type():
    analyses = [vds(5, 'gs://b/old.vds'), vds(9, 'gs://b/new.vds'), vds(12, 'gs://b/exome.vds', 'exome')]
    assert newest_vds(analyses, 'genome') == 'gs://b/new.vds'
    assert newest_vds(analyses, 'exome') == 'gs://b/exome.vds'
    assert newest_vds([], 'genome') is None


def test_vds_name_matches_cohort_demo_without_a_refresh():
    # cohort_demo's get_sg_hash(['CPG1', 'CPG2']), in any order
    assert vds_name(['CPG2', 'CPG1'], []) == '7a64ebfc1d2c25c5e17f07211230ec74d3fb13_2'


def test_refresh_changes_the_vds_name_but_not_the_count():
    plain = vds_name(['CPG1', 'CPG2'], [])
    refreshed = vds_name(['CPG1', 'CPG2'], ['CPG2'])
    assert refreshed != plain
    assert refreshed.endswith('_2')
    assert refreshed == vds_name(['CPG2', 'CPG1'], ['CPG2'])


def test_main_writes_the_sheet_and_vds_lists(monkeypatch):
    response = {
        'cohorts': [cohort('COH1', [sg('CPG2'), sg('CPG1')])],
        'project': {
            'analyses': [vds(3, 'gs://b/one.vds'), vds(4, 'gs://b/two.vds', 'exome')],
        },
    }
    queried = {}

    def fake_query(_query, variables):
        queried.update(variables)
        return response

    monkeypatch.setattr(build_samplesheet, 'query', fake_query)
    build_samplesheet.main(PROJECT, 'genome', ['COH1', 'COH1'], ['CPG2'])

    assert queried == {'project': PROJECT, 'cohorts': ['COH1']}
    with open(f'{vds_name(["CPG1", "CPG2"], ["CPG2"])}.tsv') as handle:
        assert [row['sg_id'] for row in csv.DictReader(handle, delimiter='\t')] == ['CPG1', 'CPG2']
    with open('base_vds.txt') as handle:
        assert handle.read() == 'gs://b/one.vds\n'


def test_main_writes_an_empty_base_when_there_is_none(monkeypatch):
    response = {'cohorts': [cohort('COH1', [sg('CPG1')])], 'project': {'analyses': []}}
    monkeypatch.setattr(build_samplesheet, 'query', lambda _query, variables: response)  # noqa: ARG005
    build_samplesheet.main(PROJECT, 'genome', ['COH1'], [])
    with open('base_vds.txt') as handle:
        assert handle.read() == ''


def test_every_missing_input_is_reported_at_once(missing):
    rows = build_rows([cohort('COH1', [sg('CPG1'), sg('CPG2'), sg('CPG3')])], 'genome')
    missing.update(
        {'gs://b/base.vds', 'gs://cpg-ds-main/gvcf/CPG1.g.vcf.gz', 'gs://cpg-ds-main/gvcf/CPG2.g.vcf.gz.tbi'}
    )
    with pytest.raises(FileNotFoundError) as error:
        validate_file_existence(rows, 'gs://b/base.vds')
    assert 'VDS: gs://b/base.vds does not exist' in str(error.value)
    assert 'gVCF: gs://cpg-ds-main/gvcf/CPG1.g.vcf.gz does not exist' in str(error.value)
    assert 'gVCF: gs://cpg-ds-main/gvcf/CPG2.g.vcf.gz.tbi does not exist' in str(error.value)
    assert 'CPG3' not in str(error.value)


def test_refresh_of_an_sg_not_being_combined_fails(monkeypatch):
    response = {'cohorts': [cohort('COH1', [sg('CPG1')])], 'project': {'analyses': []}}
    monkeypatch.setattr(build_samplesheet, 'query', lambda _query, variables: response)  # noqa: ARG005
    with pytest.raises(SamplesheetError, match='1 SG IDs were due to be refreshed'):
        build_samplesheet.main(PROJECT, 'genome', ['COH1'], ['CPG1', 'CPG9'])
