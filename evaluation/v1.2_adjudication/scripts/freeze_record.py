#!/usr/bin/env python3
"""Record public pre-adjudication materials; never includes the tool key."""
from __future__ import annotations
import hashlib,json,platform,subprocess,sys
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1]
FILES=['evaluation/v1.2_adjudication/frozen/benchmark_manifest.json','evaluation/v1.2_adjudication/frozen/utility_definition.md','evaluation/v1.2_adjudication/reviewer_files/blinded_adjudication_cases.csv','evaluation/v1.2_adjudication/reviewer_files/blinded_adjudication_cases.xlsx','evaluation/v1.2_adjudication/reviewer_files/adjudication_guide.md','evaluation/v1.2_adjudication/reviewer_files/adjudication_data_dictionary.md','evaluation/v1.2_adjudication/reviewer_files/error_taxonomy.md']
def sha(p): return hashlib.sha256((REPO/p).read_bytes()).hexdigest()
def main():
    data={'git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),'pr_number':15,'recorded_at_utc':datetime.now(timezone.utc).isoformat(),'benchmark_manifest_sha256':sha(FILES[0]),'utility_definition_sha256':sha(FILES[1]),'blinded_unit_file_sha256':sha(FILES[2]),'workbook_sha256':sha(FILES[3]),'reviewer_guide_sha256':sha(FILES[4]),'evidence_dictionary_sha256':sha(FILES[5]),'error_taxonomy_sha256':sha(FILES[6]),'number_of_fb_units':689,'evaluable_loci':390,'software':{'python':sys.version,'platform':platform.platform()},'confidential_materials_included':False}
    (ROOT/'frozen/pre_adjudication_freeze_record.json').write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
    print(json.dumps(data,indent=2))
if __name__=='__main__':main()
