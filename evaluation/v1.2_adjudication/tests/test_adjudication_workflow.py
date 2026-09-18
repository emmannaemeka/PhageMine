import csv, json, subprocess, sys
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
REPO=ROOT.parents[1]

def read(path):
    with path.open() as f:return list(csv.DictReader(f))

def test_review_file_is_neutral_and_complete():
    rows=read(ROOT/'reviewer_files/blinded_adjudication_cases.csv')
    assert len(rows)==689
    assert {r['case_id'] for r in rows}=={f'FB{i:04d}' for i in range(1,690)}
    assert all('PhageMine' not in r['prediction_A_evidence_raw'] and 'Pharokka' not in r['prediction_A_evidence_raw'] for r in rows)
    assert {r['review_target'] for r in rows}=={'A','B'}
    assert all(not r['final_judgment'] for r in rows)

def test_workbook_has_validation_and_neutral_sheets():
    openpyxl=pytest.importorskip('openpyxl')
    wb=openpyxl.load_workbook(ROOT/'reviewer_files/blinded_adjudication_cases.xlsx')
    assert wb.sheetnames==['Instructions','Adjudication','Decision_Definitions']
    assert len(wb['Adjudication'].data_validations.dataValidation)>=6
    assert len(list(wb['Adjudication'].rows))==690

def test_utility_definition_is_frozen():
    d=json.loads((ROOT/'frozen/utility_definition.json').read_text())
    assert d['primary_utility']['weights']['CORRECT']==1.0
    assert d['primary_utility']['weights']['PARTIALLY_CORRECT']==0.5
    assert d['primary_utility']['weights']['TOO_GENERAL']==0.25
    assert d['primary_utility']['unresolvable_rule'].startswith('missing')
    assert d['primary_margin']==0.05 and d['structural_noninferiority_margin']==-0.03 and d['unsupported_specificity_safety_bound']==0.02

def test_final_analysis_refuses_blank_or_incomplete_reviews(tmp_path):
    p=tmp_path/'blank.csv'; p.write_text('case_id,final_judgment,reviewer_id,review_timestamp\nFB0001,,,\n')
    proc=subprocess.run([sys.executable,str(ROOT/'scripts/post_adjudication.py'),'--reviews',str(p),'--output',str(tmp_path/'out'),'--final'],capture_output=True,text=True)
    assert proc.returncode!=0 and 'Reviewer provenance missing' in proc.stderr

def test_synthetic_complete_analysis_is_explicitly_separate(tmp_path):
    template=read(ROOT/'reviewer_files/blinded_adjudication_cases.csv')
    fields=list(template[0]); out=tmp_path/'synthetic.csv'
    with out.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for i,row in enumerate(template,1):
            row=dict(row);row['final_judgment']='CORRECT';row['reviewer_id']='SYNTHETIC_TEST_ONLY';row['review_timestamp']='2026-09-18T00:00:00Z';row['evidence_supports_prediction']='YES';row['evidence_strength']='STRONG';row['specificity_appropriate']='YES';row['conflicting_evidence']='NO';row['error_category']='NONE';row['reviewer_rationale']='SYNTHETIC DATA — TEST ONLY';w.writerow(row)
    result=tmp_path/'analysis'
    proc=subprocess.run([sys.executable,str(ROOT/'scripts/post_adjudication.py'),'--reviews',str(out),'--output',str(result),'--final'],capture_output=True,text=True)
    assert proc.returncode==0
    data=json.loads((result/'statistical_report.json').read_text()); assert data['status']=='ADJUDICATED'
    assert data['definition']['primary_margin']==0.05
    assert data['genome_results']
    assert not (ROOT/'analysis'/'statistical_report.json').exists()

def test_error_analysis_requires_unblind_flag(tmp_path):
    p=tmp_path/'review.csv'; p.write_text('case_id,final_judgment,reviewer_id,review_timestamp\n')
    proc=subprocess.run([sys.executable,str(ROOT/'scripts/error_analysis.py'),'--reviews',str(p),'--output',str(tmp_path/'out')],capture_output=True,text=True)
    assert proc.returncode!=0 and '--unblind' in proc.stderr

def test_frozen_manifest_and_checksums():
    manifest=json.loads((ROOT/'frozen/benchmark_manifest.json').read_text())
    assert manifest['repository_commit']=='4933d8fa532a0157c24c3c01dbf1d8adfa787ed7'
    assert manifest['counts']['PhageMine']['named_assertions']==358
    assert manifest['counts']['Pharokka']['unresolved']==38
    assert manifest['named_yield']==[144,56,50,95,8,48,38]
