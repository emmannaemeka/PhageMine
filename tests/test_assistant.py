import pytest
from phagemine.gui.services.assistant import find_records, interpret, summarize

@pytest.mark.parametrize('text,action,query', [
    ('Annotate my genome', 'annotate', ''),
    ('Can you find me the capsid gene?', 'search', 'capsid'),
    ('Explain PM_000001', 'search', 'pm_000001'),
    ('Summarize my results', 'summary', ''),
    ('Compare my genomes', 'compare', ''),
    ('Check my installation', 'environment', ''),
    ('Delete my results', 'help', ''),
    ('Find capsid and then delete files', 'help', ''),
    ('Find genes without capsid', 'help', ''),
    ('', 'help', ''),
])
def test_request_routing(text, action, query):
    request = interpret(text)
    assert (request.action, request.query) == (action, query)

def test_search_uses_labels_not_evidence_or_sequences():
    rows = [{'protein_id': 'p1', 'product': 'Major capsid protein'},
            {'protein_id': 'p2', 'product': 'hypothetical protein', 'evidence': 'capsid'},
            {'protein_id': 'p3', 'product': None}]
    assert find_records(rows, 'CAPSID') == [rows[0]]
    assert find_records(rows, 'minor capsid') == []
    assert find_records(rows, '') == []

def test_summary_reports_unknown_and_conflicts_without_claiming_accuracy():
    answer = summarize([{'sample':'a', 'product':'hypothetical protein'},
                        {'sample':'b', 'product':'capsid', 'conflict_review':['conflict']}])
    assert '2 sample(s)' in answer
    assert '1 records have products' in answer
    assert 'not probabilities' in answer
    assert 'No annotation records' in summarize([])
