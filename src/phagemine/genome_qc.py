"""Conservative sequence QC and Turner et al. (2021) review checklist.

Findings are candidates for review, never proof of physical termini or safety.
Coordinates are 1-based inclusive. This module does not alter sequence or IDs.
"""
from __future__ import annotations
import hashlib
import json
import re
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

GUIDE = 'https://pmc.ncbi.nlm.nih.gov/articles/PMC8785237/'
QUESTIONS = [
    'Sequencing and coverage', 'Genome completeness and physical termini',
    'Ambiguous bases and potential frameshifts', 'Co-linearity and assembly artefacts',
    'Structural annotation and small ORFs', 'CDS start-site curation',
    'Functional evidence', 'Lysogeny, toxins and AMR screening',
    'Regulatory elements and noncoding RNA', 'Isolate/prophage/environmental source',
    'Phage name and locus tags', 'Submission cleanup', 'Taxonomic assignment',
    'Database updates and reproducibility',
]


def sanitize_product(value, supported=True):
    value = ' '.join(str(value or '').split())
    if not supported or not value:
        return 'hypothetical protein'
    vague = re.search(r'\b(hypothetical|uncharacteri[sz]ed|unknown|unnamed)\b', value, re.I)
    identifier = re.search(r'\b(?:[A-Z]{1,4}_?\d{4,}(?:\.\d+)?|(?:ORF|gp|PHROG|VOG|PF)\d+)\b', value, re.I)
    symbol = re.fullmatch(r'(?:[a-z]{2,3}[A-Z][0-9]?|[A-Z][a-z]{2}[A-Z]?)', value)
    if vague or identifier or symbol or value.lower() in {'putative phage protein', 'phage protein', 'putative protein'}:
        return 'hypothetical protein'
    return value


def overlap(a, b):
    return max(0, min(a['end'], b['end']) - max(a['start'], b['start']) + 1)


def sequence_findings(sequence, cds, homopolymer_min=10):
    findings = []
    patterns = [('AMBIGUOUS_BASES', r'[^ACGT]+'), ('HOMOPOLYMER', rf'(A{{{homopolymer_min},}}|C{{{homopolymer_min},}}|G{{{homopolymer_min},}}|T{{{homopolymer_min},}})')]
    for code, pattern in patterns:
        for m in re.finditer(pattern, sequence.upper()):
            row = dict(code=code, start=m.start()+1, end=m.end(), sequence=m.group())
            row['cds_ids'] = [p['id'] for p in cds if overlap(row, p)]
            row['invalid_symbols'] = sorted(set(m.group()) - set('ACGTRYSWKMBDHVN'))
            findings.append(row)
    return findings


def terminal_repeats(sequence, minimum=20, maximum=10000):
    """Exact end repeats, bounded to half the contig; no topology inference."""
    seq = sequence.upper(); results = []
    for kind in ('DTR', 'ITR'):
        for size in range(min(maximum, len(seq)//2), minimum-1, -1):
            left, right = seq[:size], seq[-size:]
            if set(left+right) - set('ACGT'): continue
            if left == (right if kind == 'DTR' else str(Seq(right).reverse_complement())):
                internal = []; pos = seq.find(left, size)
                while pos >= 0 and pos + size <= len(seq)-size:
                    internal.append(pos+1); pos=seq.find(left,pos+1)
                results.append({'code': kind+'_CANDIDATE', 'length_nt': size, 'internal_copies': internal})
                break
    return results


def frameshift_candidates(cds, hits, maximum_gap=100):
    """Adjacent CDSs sharing a subject, different frames, ordered subject pieces.

    hits: id -> [{target, subject_start, subject_end, source, evalue}]. Domain
    identifiers without aligned subject coordinates cannot establish a split hit.
    """
    findings=[]
    for strand in ('+', '-'):
        genes=sorted((p for p in cds if p['strand']==strand), key=lambda p:p['start'], reverse=strand=='-')
        for a,b in zip(genes,genes[1:]):
            gap=(b['start']-a['end']-1) if strand=='+' else (a['start']-b['end']-1)
            frame_a=(a['start']-1)%3 if strand=='+' else a['end']%3
            frame_b=(b['start']-1)%3 if strand=='+' else b['end']%3
            if not -30<=gap<=maximum_gap or frame_a==frame_b: continue
            for x in hits.get(a['id'],[]):
                for y in hits.get(b['id'],[]):
                    if x.get('target')!=y.get('target') or x.get('source')!=y.get('source'): continue
                    if max(x.get('evalue',1),y.get('evalue',1))>1e-5: continue
                    if not all(k in h for h in (x,y) for k in ('subject_start','subject_end')): continue
                    if x['subject_start'] < y['subject_start'] and x['subject_end'] <= y['subject_start']+10:
                        findings.append({'code':'PUTATIVE_FRAMESHIFT', 'cds_ids':[a['id'],b['id']], 'target':x['target'], 'gap_nt':gap, 'reason':'Ordered split homology in different frames; sequencing error, programmed frameshift or intron requires inspection'})
    return findings


def safety_markers(proteins):
    patterns={
        'lysogeny':r'\b(integrase|excisionase|site.specific recombinase|CI.like repressor|lysogen\w*)\b',
        'toxin_virulence':r'\b(toxin|enterotoxin|exotoxin|cytotoxin|shiga|virulence)\b',
        'amr':r'\b(beta.lactamase|antimicrobial resistance|antibiotic resistance|aminoglycoside.*(?:acetyl|phospho|adenylyl)transferase|tetracycline resistance)\b',
    }
    findings=[]
    for p in proteins:
        for e in p.get('evidence',[]):
            if str(e.get('status','')).upper()!='REAL' or not e.get('supports',False): continue
            text=' '.join(str(e.get(k) or '') for k in ('statement','description','family_name'))
            for category,pattern in patterns.items():
                if re.search(pattern,text,re.I):
                    # Anti-toxins are not themselves evidence for a toxin gene.
                    if category=='toxin_virulence' and re.search(r'\banti[- ]?toxin\b',text,re.I) and not re.search(r'\b(?:shiga|entero|exo|cyto)toxin\b',text,re.I): continue
                    findings.append({'category':category,'id':p['protein_id'],'source':e.get('source'),'identifier':e.get('identifier'),'evidence':text})
    return {'status':'REVIEW' if findings else 'NO_MARKERS_IN_AVAILABLE_EVIDENCE', 'findings':findings,
            'dedicated_screening':'NOT_ASSESSED', 'therapeutic_safety':'NOT_ESTABLISHED',
            'limitation':'Evidence-text screening is not a dedicated AMR/VFDB search; absence of hits is not evidence of safety.'}


def curate_regulatory(features, cds):
    accepted=[]; rejected=[]
    for f in features:
        if any(overlap(f,p) for p in cds):
            rejected.append({**f,'reason':'REGULATORY_FEATURE_OVERLAPS_CDS'})
        else:
            accepted.append({**f,'inference':'ab initio prediction','note':'putative; computational prediction'})
    return accepted,rejected


def write_qc(output, genome=None, proteins=None, metadata=None, provenance=None):
    """Inspect final exported GenBank, including any subsequently added RNA."""
    root=Path(output); package=root/'genbank_submission'
    if not package.exists(): package=root
    gb=package/'genome.gbk'
    record=SeqIO.read(gb,'genbank')
    sequence=str(record.seq).upper(); cds=[]; rna=[]; regulatory=[]
    for f in record.features:
        row={'id':f.qualifiers.get('locus_tag',['unknown'])[0], 'start':int(f.location.start)+1,'end':int(f.location.end),'strand':'+' if f.location.strand==1 else '-'}
        if f.type=='CDS':
            nt=str(f.extract(record.seq)).upper()
            partial5=type(f.location.start if f.location.strand==1 else f.location.end).__name__!='ExactPosition'
            partial3=type(f.location.end if f.location.strand==1 else f.location.start).__name__!='ExactPosition'
            offset=int(f.qualifiers.get('codon_start',['1'])[0])-1
            coding=nt[offset:]; translated=str(Seq(coding[:len(coding)//3*3]).translate(table=11))
            from .start_selection import STARTS_11, STOPS
            errors=[]
            if not partial5 and coding[:3] not in STARTS_11: errors.append('CDS_START_INVALID')
            if not partial3 and coding[-3:] not in STOPS: errors.append('CDS_STOP_MISSING')
            if len(coding)%3 and not partial3: errors.append('CDS_FRAME_INVALID')
            if '*' in translated[:-1]: errors.append('CDS_INTERNAL_STOP')
            if f.qualifiers.get('transl_table')!=['11']: errors.append('TRANSL_TABLE_NOT_11')
            row.update(errors=errors,start_codon=coding[:3],partial_5prime=partial5,partial_3prime=partial3,product=f.qualifiers.get('product',[''])[0])
            cds.append(row)
        elif f.type in {'tRNA','tmRNA','rRNA','ncRNA'}: rna.append(row)
        elif f.type in {'promoter','terminator'}: regulatory.append(row)
    def read(path, default):
        return json.loads(path.read_text()) if path.exists() else default
    if proteins is None:
        evidence=read(root/'evidence.json',{})
        proteins=evidence.get('proteins',[]) if isinstance(evidence,dict) else evidence
    proteins=[p.to_dict() if hasattr(p,'to_dict') else p for p in proteins]
    metadata=metadata or read(package/'submission_metadata.json',{})
    provenance=provenance or read(package/'provenance.json',{})
    if hasattr(metadata,'__dataclass_fields__'):
        from dataclasses import asdict
        metadata=asdict(metadata)
    anomalies=sequence_findings(sequence,cds); repeats=terminal_repeats(sequence)
    hits={}
    for p in proteins:
        for e in p.get('evidence',[]):
            m=dict(e.get('metrics',{}))
            if 'subject_coordinates' in m:
                m.update(subject_start=m['subject_coordinates']['start'], subject_end=m['subject_coordinates']['end'])
            if str(e.get('status','')).upper()=='REAL' and 'subject_start' in m and 'subject_end' in m:
                hits.setdefault(p['protein_id'],[]).append({**m,'target':e.get('identifier'),'source':e.get('source')})
    shifts=frameshift_candidates(cds,hits)
    safety=safety_markers(proteins)
    dedicated=read(root/'safety_screen.json',{})
    if dedicated:
        query=package/'proteins.faa'
        if query.exists() and dedicated.get('query_sha256')==hashlib.sha256(query.read_bytes()).hexdigest():
            safety['dedicated_screening']=dedicated
            safety['findings'].extend(dedicated.get('findings',[]))
        else:
            safety['dedicated_screening']='STALE_QUERY_MISMATCH'
    if safety['findings']: safety['status']='REVIEW'

    overlaps=[{'rna_id':a['id'],'cds_id':b['id'],'overlap_nt':overlap(a,b)} for a in rna for b in cds if overlap(a,b)]
    _,bad_reg=curate_regulatory(regulatory,cds)
    tags=[p['id'] for p in cds]
    bad_tags=[t for t in tags if not re.fullmatch(r'[A-Za-z][A-Za-z0-9]*_\d+',t)]
    if tags and not bad_tags:
        numbers=sorted(int(t.rsplit('_',1)[1]) for t in tags)
        if numbers != list(range(numbers[0],numbers[0]+len(numbers))):
            bad_tags.append('NONSEQUENTIAL_LOCUS_NUMBERS')
    sources={e.get('source') for p in proteins for e in p.get('evidence',[]) if e.get('supports') and str(e.get('status','')).upper()=='REAL'}
    criteria=[]
    def item(n,status,reason,details=None):
        criteria.append({'question':n,'criterion':QUESTIONS[n-1],'status':status,'reason':reason,'details':details or []})
    item(1,'WARN','Read mapping, coverage uniformity and library preparation require supplied read evidence')
    item(2,'WARN','Physical completeness is not established by sequence inspection',repeats)
    circular=metadata.get('topology')=='circular' or record.annotations.get('topology')=='circular'
    item(3,'WARN' if anomalies or shifts else 'PASS','Sequence anomaly and available split-hit screen; missing alignments cannot exclude frameshifts',anomalies+shifts)
    item(4,'WARN','Circular topology requires physical verification' if circular else 'Co-linearity requires a related genome comparison',repeats)
    refinement=read(root/'gene_calls/start_refinement/manifest.json',{})
    item(5,'PASS' if refinement else 'WARN','Intergenic ORF evaluation recorded' if refinement else 'No intergenic rescue audit available',[refinement] if refinement else [])
    bad=[p for p in cds if p['errors']]
    partial=[p for p in cds if p['partial_5prime'] or p['partial_3prime']]
    item(6,'FAIL' if bad else 'WARN' if partial else 'PASS','Codon/frame checks; biological start positions remain computational',bad+partial)
    evidence_audit=[]
    for protein in proteins:
        real=[e for e in protein.get('evidence',[]) if e.get('supports') and str(e.get('status','')).upper()=='REAL']
        modalities={e.get('modality') for e in real}
        evidence_audit.append({'id':protein['protein_id'],'sources':sorted({e.get('source','unknown') for e in real}),
            'sequence_and_profile':bool(modalities & {'sequence_similarity','phage_orthology'}) and bool(modalities & {'domain','phage_profile_hmm','viral_orthology'})})
    item(7,'PASS' if evidence_audit and all(e['sequence_and_profile'] for e in evidence_audit) else 'WARN',
         'Per-protein evidence diversity; weak functional assignments remain hypothetical',evidence_audit)
    item(8,'WARN','Therapeutic suitability is not established by these checks',[safety])
    item(9,'FAIL' if bad_reg else 'WARN','RNA overlaps need review; tmRNA search and regulatory predictions may be unconfigured',overlaps+bad_reg)
    missing=[k for k in ('organism','host','isolation_source','authors','submitter_name','submitter_email') if not metadata.get(k)]
    item(10,'WARN' if missing else 'PASS','Source metadata completeness',missing)
    item(11,'FAIL' if bad_tags or len(set(tags))!=len(tags) else 'WARN','Local locus syntax checked; prefix registration and name uniqueness require confirmation',bad_tags)
    bad_products=[p for p in cds if sanitize_product(p['product'])!=p['product']]
    item(12,'FAIL' if bad_products else 'PASS','Exported product-name cleanup',bad_products)
    item(13,'WARN','Taxonomic placement requires independent review against current taxonomy')
    item(14,'PASS' if provenance else 'WARN','Provenance retained; the submitter remains responsible for later updates')
    status='FAIL' if any(c['status']=='FAIL' for c in criteria) else 'WARN' if any(c['status']=='WARN' for c in criteria) else 'PASS'
    report={'schema_version':1,'guide':GUIDE,'status':status,'insdc_acceptance':False,'sequence_sha256':hashlib.sha256(sequence.encode()).hexdigest(),'coordinate_system':'1-based-inclusive','criteria':criteria,'safety':safety,'cds_count':len(cds),'rna_count':len(rna),'provenance':provenance}
    for dest in {root,package}:
        (dest/'phagemine_qc_report.json').write_text(json.dumps(report,indent=2,default=str)+'\n')
        lines=['# PhageMine genome QC',f'Overall: **{status}**. This is not INSDC acceptance or therapeutic clearance.',f'Guide: {GUIDE}','','| Question | Check | Status | Finding |','|---|---|---|---|']
        lines += [f"| {c['question']} | {c['criterion']} | {c['status']} | {c['reason']} |" for c in criteria]
        lines += ['', '## Findings', '```json', json.dumps({str(c['question']):c['details'] for c in criteria if c['details']},indent=2,default=str), '```']
        (dest/'phagemine_qc_report.md').write_text('\n'.join(lines)+'\n')
    return report
