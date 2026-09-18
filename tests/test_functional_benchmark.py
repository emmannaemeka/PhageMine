import csv
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import pytest
from phagemine.functional_benchmark import (CATEGORIES, PENDING, classify, informative,
    load_synonyms, match_loci, paired_test, summarize, wilson)

ROOT=Path(__file__).resolve().parents[1]
PACKAGE=ROOT/'evaluation/v1.2_functional'


def module(path):
    spec=importlib.util.spec_from_file_location(path.stem,path)
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


def read(path):
    with path.open() as f: return list(csv.DictReader(f,delimiter='\t'))


def test_published_aggregate_consistency():
    plot=module(ROOT/'docs/benchmark_v1.2/scripts/make_figures.py')
    totals=plot.load_totals(ROOT/'docs/benchmark_v1.2')
    per_genome=read(ROOT/'docs/benchmark_v1.2/tables/functional_yield_by_genome.tsv')
    for row in totals:
        tool=row['tool']; n=int(row['named_product_calls'])
        assert n=={'PhageMine':439,'Pharokka':398,'Prokka':328}[tool]
        assert sum(int(r['named_product_calls']) for r in per_genome if r['tool']==tool)==n
        for doc in ['README.md','docs/BENCHMARK.md']:
            line=next(l for l in (ROOT/doc).read_text().splitlines() if l.startswith('| '+tool+' |'))
            cells=[x.strip() for x in line.split('|')[1:-1]]
            assert int(cells[1])==int(row['predicted_cds']) and int(cells[-1])==n
            assert cells[2:4]==(['0.8846','0.9352'] if tool=='Prokka' else ['0.7762','0.8821'])


def test_plot_independent_of_invalid_per_genome_zeros(tmp_path):
    plot=module(ROOT/'docs/benchmark_v1.2/scripts/make_figures.py')
    (tmp_path/'tables').mkdir()
    (tmp_path/'tables/functional_yield_totals.tsv').write_bytes((ROOT/'docs/benchmark_v1.2/tables/functional_yield_totals.tsv').read_bytes())
    (tmp_path/'tables/functional_yield_by_genome.tsv').write_text('tool\tnamed_product_calls\nPhageMine\t0\n')
    assert int(plot.load_totals(tmp_path)[0]['named_product_calls'])==439
    p=tmp_path/'tables/functional_yield_totals.tsv'; p.write_text(p.read_text().replace('439','0'))
    with pytest.raises(AssertionError): plot.load_totals(tmp_path)


@pytest.mark.parametrize('text',['hypothetical protein','hypothetical protein containing ASCH domain',
    'Uncharacterised protein','uncharacterized 8.1 kDa protein','unknown function',
    'conserved protein of unknown function','conserved phage protein','protein',''])
def test_uninformative_labels(text):
    assert not informative(text)
    assert classify('endolysin',text,{})=='NON_EVALUABLE_REFERENCE'
    assert classify(text,'endolysin',{})=='ABSTENTION'


def test_normalization_and_explicit_synonyms():
    rules=load_synonyms(PACKAGE/'synonym_rules.tsv')
    assert classify('Major capsid protein','major capsid protein',rules)=='EXACT_PRODUCT_AGREEMENT'
    assert classify('major head protein','major capsid protein',rules)=='EQUIVALENT_FUNCTION'
    assert classify('putative major capsid protein','major capsid protein',rules)==PENDING
    assert classify('DNA polymerase','DNA helicase',rules)==PENDING
    assert classify('terminase small subunit','terminase large subunit',rules)==PENDING
    assert classify('major head protein','major capsid protein',{})==PENDING


def test_bad_synonyms_rejected(tmp_path):
    p=tmp_path/'rules.tsv'; p.write_text('rule_id\tterm_a\tterm_b\trationale\nx\thypothetical protein\tcapsid protein\ttest\n')
    with pytest.raises(ValueError): load_synonyms(p)


def rec(start=1,end=100,accession='a',strand='+'):
    return dict(start=start,end=end,accession=accession,strand=strand)


def test_coordinate_accession_and_strand_matching():
    assert match_loci([rec()],[rec()])[0]=={0:0}
    for ref in [rec(accession='b'),rec(strand='-'),rec(end=101)]:
        assert match_loci([rec()],[ref])[0]=={}
    assert match_loci([rec()],[rec(end=101)],'relaxed')[0]=={0:0}
    with pytest.raises(ValueError): match_loci([rec(),rec()],[rec()])


def test_ambiguous_overlap_never_forced():
    matches,ambiguous=match_loci([rec()],[rec(2,90),rec(3,95)],'relaxed')
    assert not matches and ambiguous=={0:[0,1]}
    matches,ambiguous=match_loci([rec(),rec(2,95)],[rec(3,90)],'relaxed')
    assert not matches and len(ambiguous)==2
    assert match_loci([rec(),rec(2,95)],[rec()],'relaxed')[0]=={0:0}


def test_categories_metrics_and_non_evaluable_protection():
    result=summarize(list(CATEGORIES))
    assert result['evaluable_named_reference_loci']==7
    assert result['named_assertions']==6
    assert result['strict_functional_precision']==2/6
    assert result['expanded_compatible_precision']==4/6
    assert result['functional_recall']==2/7
    assert result['unsupported_specificity_rate']==1/6
    assert summarize(['NON_EVALUABLE_REFERENCE'])['strict_functional_precision'] is None
    assert classify('hypothetical protein','hypothetical protein',{})=='NON_EVALUABLE_REFERENCE'


def test_extra_named_call_is_not_automatically_correct():
    result=summarize([classify('endolysin','DNA helicase',{})])
    assert result['unresolved']==1
    assert result['strict_functional_precision'] is None
    assert result['functional_coverage']==1
    assert result['strict_functional_precision_identification_low']==0
    assert result['strict_functional_precision_identification_high']==1


def test_paired_statistics_and_intervals():
    assert paired_test([PENDING],['ABSTENTION'])['p_value'] is None
    assert paired_test(['EXACT_PRODUCT_AGREEMENT']*6,['ABSTENTION']*6)['p_value']==.03125
    assert wilson(0,0)==(None,None)
    assert wilson(0,10)[0]==0


def test_blinding_reproducible_and_reversible():
    run=module(PACKAGE/'scripts/run_benchmark.py')
    rows=read(PACKAGE/'functional_by_locus.tsv')
    a,key=run.blinded_rows(rows)
    assert (a,key)==run.blinded_rows(rows)
    assert {k['prediction_A_tool'] for k in key}=={'PhageMine','Pharokka'}
    assert all('tool' not in field for field in a[0])
    assert len({r['adjudication_id'] for r in a})==len(a)
    indexed={(r['matching_mode'],r['accession'],int(r['start']),int(r['end']),r['strand']):r for r in rows}
    for blinded,mapping in zip(a,key):
        r=indexed[(blinded['matching_mode'],blinded['accession'],int(blinded['start']),int(blinded['end']),blinded['strand'])]
        assert blinded['prediction_A']==r[mapping['prediction_A_tool']+'_product']
        assert blinded['prediction_B']==r[mapping['prediction_B_tool']+'_product']


def test_review_ingestion_and_tamper_protection(tmp_path):
    run=module(PACKAGE/'scripts/run_benchmark.py')
    row=dict(rec(),matching_mode='exact',reference_start=1,reference_end=100,reference_product='endolysin',
        PhageMine_product='lysin',Pharokka_product='endolysin',PhageMine_category=PENDING,
        Pharokka_category='EXACT_PRODUCT_AGREEMENT',PhageMine_review_notes='',Pharokka_review_notes='')
    review,key=run.blinded_rows([row]); submitted=dict(review[0])
    side='A' if key[0]['prediction_A_tool']=='PhageMine' else 'B'
    submitted['reviewer_category_'+side]='COMPATIBLE_BUT_MORE_SPECIFIC'; submitted['reviewer_notes']='asserted'
    path=tmp_path/'review.tsv';run.write(path,[submitted])
    with pytest.raises(ValueError,match='evidence'):run.apply_reviews([row],review,key,path)
    submitted['reviewer_category_'+side]='EQUIVALENT_FUNCTION';run.write(path,[submitted])
    run.apply_reviews([row],review,key,path)
    assert row['PhageMine_category']=='EQUIVALENT_FUNCTION'
    submitted['prediction_A']='tampered';run.write(path,[submitted])
    with pytest.raises(ValueError,match='changed'):run.apply_reviews([row],review,key,path)


def test_frozen_reproduction(tmp_path):
    pytest.importorskip('Bio')
    subprocess.run([sys.executable,str(PACKAGE/'scripts/run_benchmark.py'),'--output',str(tmp_path)],check=True,capture_output=True)
    for filename in ['functional_benchmark.json','functional_summary.tsv','functional_truth_by_locus.tsv',
                     'manual_adjudication_blinded.tsv','manual_adjudication_key.tsv',
                     'phagemine_only_functional_assignments.tsv','pharokka_only_functional_assignments.tsv']:
        assert (tmp_path/filename).read_bytes()==(PACKAGE/filename).read_bytes()
    result=json.loads((tmp_path/'functional_benchmark.json').read_text())
    assert result['exclusive_assignments']['PhageMine']['total']==75
    assert result['exclusive_assignments']['Pharokka']['total']==34
    assert result['summaries'][0]['matched_reference_cds']==588


def test_actual_plot_reads_aggregate_with_invalid_per_genome(tmp_path, monkeypatch):
    pytest.importorskip('matplotlib')
    import matplotlib.pyplot as plt
    plot=module(ROOT/'docs/benchmark_v1.2/scripts/make_figures.py')
    (tmp_path/'tables').mkdir();(tmp_path/'figures').mkdir()
    (tmp_path/'tables/functional_yield_totals.tsv').write_bytes((ROOT/'docs/benchmark_v1.2/tables/functional_yield_totals.tsv').read_bytes())
    (tmp_path/'tables/functional_yield_by_genome.tsv').write_text('tool\tnamed_product_calls\nPhageMine\t0\n')
    plot.__file__=str(tmp_path/'scripts/make_figures.py')
    original=plt.Axes.bar; observed=[]
    def capture(self,x,height,*args,**kwargs):
        observed.extend(height)
        return original(self,x,height,*args,**kwargs)
    monkeypatch.setattr(plt.Axes,'bar',capture)
    plot.main()
    assert observed==[439,398,328]
    assert all((tmp_path/f'figures/functional_yield.{ext}').stat().st_size>0 for ext in ['png','svg','pdf'])


def test_archived_evidence_and_package_checksums():
    import hashlib
    for folder in [PACKAGE,ROOT/'docs/benchmark_v1.2']:
        for line in (folder/'SHA256SUMS').read_text().splitlines():
            checksum,path=line.split('  ',1)
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==checksum,path
    for row in read(PACKAGE/'source_provenance.tsv'):
        assert hashlib.sha256((ROOT/row['archived_path']).read_bytes()).hexdigest()==row['sha256']
    truth=read(PACKAGE/'functional_truth_by_locus.tsv')
    assert sum(r['compound_location']=='True' for r in truth)==11
    assert len(truth)==711
