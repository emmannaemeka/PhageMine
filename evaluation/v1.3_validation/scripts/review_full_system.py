#!/usr/bin/env python3
import sys,json,csv,hashlib,collections,shutil,argparse
from pathlib import Path
root=Path(__file__).resolve().parents[3];sys.path.insert(0,str(root/'src'))
from phagemine.resume import reclassify
from phagemine.benchmark import import_phagemine
from phagemine.curated_benchmark import evaluate
from phagemine.operational_validation import validate_run
from run_full_system import verify_full_run, REQUIRED_RESOURCES
parser=argparse.ArgumentParser(description='Audit current decision rules against fresh full-system evidence without rerunning searches')
parser.add_argument('--source',required=True,type=Path)
parser.add_argument('--output',required=True,type=Path)
args=parser.parse_args();raw=args.source.resolve();out=args.output.resolve()
if out.exists(): parser.error('Output must be a new directory')
manifest=json.loads((raw/'manifest.json').read_text())
if manifest.get('fresh_annotation_run') is not True or len(manifest['genomes']) != 10: raise ValueError('Expected the fresh ten-genome full-system run')
if set(manifest.get('required_resources', [])) != REQUIRED_RESOURCES: raise ValueError('Expected all six resources')
for row in csv.DictReader((raw/'SHA256SUMS.tsv').open(), delimiter='\t'):
 source=(raw/row['path']).resolve()
 if not source.is_relative_to(raw): raise ValueError('Checksum path escapes source directory')
 if hashlib.sha256(source.read_bytes()).hexdigest() != row['sha256']: raise ValueError('Source checksum mismatch: '+row['path'])
out.mkdir(parents=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def write(p,rows):
 with p.open('w',newline='')as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
rows=[];audits=[];changes=[]
for g in json.loads((raw/'manifest.json').read_text())['genomes']:
 accession=g['accession'];src=raw/'runs'/accession;dest=out/'runs'/accession
 verify_full_run(src)
 reclassify(src,dest)
 integrity=validate_run(dest)
 before=import_phagemine(src);after=import_phagemine(dest)
 if [(r['start'],r['end'],r['strand'])for r in before] != [(r['start'],r['end'],r['strand'])for r in after]: raise ValueError('Gene coordinates changed')
 if json.loads((src/'evidence.json').read_text()) != json.loads((dest/'evidence.json').read_text()): raise ValueError('Evidence objects changed')
 audits.append({'accession':accession,'raw_evidence_sha256':sha(src/'evidence.json'),'updated_evidence_sha256':sha(dest/'evidence.json'),'evidence_objects_unchanged':True,'coordinates_unchanged':True,'integrity':integrity})
 for a,b in zip(before,after):
  row={'accession':accession,'start':b['start'],'end':b['end'],'strand':b['strand'],'product':b.get('product')or'hypothetical protein','gene':b.get('gene')or''};rows.append(row)
  if a.get('product')!=b.get('product') or a.get('gene')!=b.get('gene'):changes.append({**row,'raw_product':a.get('product')or'', 'raw_gene':a.get('gene')or''})
 print(accession,'validated',flush=True)
write(out/'PhageMine-full-current.tsv',rows)
if changes: write(out/'decision_changes.tsv',changes)
else: (out/'decision_changes.tsv').write_text('accession\tstart\tend\tstrand\tproduct\tgene\traw_product\traw_gene\n')
(out/'reclassification_audit.json').write_text(json.dumps(audits,indent=2)+'\n')
for panel in ['original','external','combined']:
 d=out/'evaluation'/panel;d.mkdir(parents=True,exist_ok=True)
 accessions={g['accession']for g in json.loads((raw/'manifest.json').read_text())['genomes']if panel=='combined'or g['panel']==panel}
 write(d/'predictions.tsv',[r for r in rows if r['accession']in accessions]);shutil.copyfile(raw/f'evaluation/{panel}/references.tsv',d/'references.tsv')
 s=evaluate(d/'predictions.tsv',d/'references.tsv',d/'scores',synonyms=root/'evaluation/v1.3_validation/naming_equivalences.tsv')['summary']
 print(panel,{k:s[k]for k in ['predicted_models','exact_model_matches','named_assertions','EXACT_PRODUCT_AGREEMENT','EQUIVALENT_FUNCTION','unresolved','abstentions_on_exact_named_reference_loci']},flush=True)
print('changed',len(changes))
