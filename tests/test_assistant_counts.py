"""Regression tests for unsupported small counts and cross-sample count swaps."""
import json

import pytest

from phagemine.assistant import agent, grounding
from phagemine.assistant.llm import ScriptedProvider
from phagemine.assistant.store import Project


@pytest.fixture
def project(tmp_path, monkeypatch):
    monkeypatch.setenv('PHAGEMINE_PROJECTS', str(tmp_path / 'projects'))
    return Project.create('count checks')


@pytest.mark.parametrize('count', ['0', '1', '7', '9', 'seven'])
def test_unsupported_small_counts(project, count):
    assert grounding.check(f'There are {count} shared genes.', '{"shared_genes": 2}', project)


@pytest.mark.parametrize('count', [0, 1, 7, 9])
def test_supported_small_counts(project, count):
    assert grounding.check(f'There are {count} shared genes.', json.dumps({'shared_genes': count}), project) == []


@pytest.mark.parametrize('corpus', [
    '{"genome_A_genes":10,"genome_B_genes":20}',
    '{"sample":"Genome A","proteins":10}\n{"sample":"Genome B","proteins":20}',
    '{"proteins_with_homolog_in_other":{"Genome A":10,"Genome B":20}}',
])
@pytest.mark.parametrize('answer', [
    'Genome A has 20 genes and Genome B has 10 genes.',
    '| Genome A | 20 |\n| Genome B | 10 |',
    'Genome A: 20, Genome B: 10.',
])
def test_swapped_sample_counts_are_rejected(project, corpus, answer):
    assert any('sample' in p for p in grounding.check(answer, corpus, project))


def test_correct_sample_counts_are_accepted(project):
    corpus = '{"proteins_with_homolog_in_other":{"phageA":2,"phageB":7}}'
    assert grounding.check('phageA: 2, phageB: 7.', corpus, project) == []
    assert grounding.check('| phageA | 2 |\n| phageB | 7 |', corpus, project) == []


def test_other_numeric_fields_do_not_authorise_another_sample(project):
    corpus = '{"sample":"phageA","proteins":2}\n{"sample":"phageB","proteins":7,"score":99}'
    assert grounding.check('phageA has 7 proteins.', corpus, project)


def test_presentation_numbers_are_not_counts(project):
    assert grounding.check('1. Inspect the results.\n2. Continue to step 3.', '{}', project) == []
    assert grounding.check('1. There are 7 shared genes.', '{"shared_genes":2}', project)


def test_sample_boundaries_do_not_match_prefixes(project):
    corpus = '{"proteins_with_homolog_in_other":{"phage":2,"phage-long":7}}'
    assert grounding.check('phage-long: 7.', corpus, project) == []
    assert grounding.check('phage: 7.', corpus, project)


def test_error_result_does_not_authorise_swapped_count(project):
    corpus = '{"sample":"phageA","proteins":2}\n{"ok":false,"sample":"phageA","proteins":7}'
    assert grounding.check('phageA has 7 proteins.', corpus, project)


def test_retry_withholds_unsupported_small_count(project):
    provider = ScriptedProvider([[{'type':'text', 'text':'There are 7 shared genes.'}]] * 2)
    result = agent.run_turn(project, 'How many genes are shared?', provider)
    assert len(provider.calls) == 2
    assert result['grounding_problems']
    assert '7 shared genes' not in result['text']
    assert 'withheld' in result['text']


def test_number_inside_sample_label_is_not_a_claim(project):
    corpus = '{"proteins_with_homolog_in_other":{"genome_7":2}}'
    assert grounding.check('genome 7 has 2 proteins.', corpus, project) == []


def test_agent_withholds_cross_sample_swap(project, monkeypatch):
    monkeypatch.setattr(agent.tools, 'execute', lambda *a: {
        'ok': True, 'proteins_with_homolog_in_other': {'phageA': 2, 'phageB': 7}})
    answer = [{'type': 'text', 'text': 'phageA has 7 proteins and phageB has 2 proteins.'}]
    provider = ScriptedProvider([
        [{'type': 'tool_use', 'id': 'count_lookup', 'name': 'get_shared_genes', 'input': {}}],
        answer, answer,
    ])
    result = agent.run_turn(project, 'Compare the counts.', provider)
    assert result['grounding_problems']
    assert 'phageA has 7' not in result['text']
    assert 'withheld' in result['text']
