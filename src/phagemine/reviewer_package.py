"""Self-contained final exports and a lossless, filterable reviewer workbook."""
from pathlib import Path
import json
import shutil
import collections
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import BeforePosition,AfterPosition
from .start_refinement import rbs_score

COLUMNS=['Locus_Tag','Start','End','Strand','Length_nt','Length_aa','Partial_5_Prime','Partial_3_Prime','Codon_Start','Start_Codon','Stop_Codon','Rare_Start_Flag','Upstream_20bp_DNA','RBS_Motif','RBS_Spacer_nt','RBS_Score','Start_Total_Score','Start_Selection_Reason','Caller_Support','Confidence_Tier','Adjacent_Gene_Overlap_nt','Product_Name','Top_PHROGs_Hit','PHROGs_Evalue','Top_VOGDB_Hit','VOGDB_Evalue','Top_SwissProt_Hit','SwissProt_Evalue','Domain_Hits','Codon_Verifier_Status','Review_Flags','Protein_Translation','Nucleotide_Sequence']


def read(path, default):
    return json.loads(path.read_text()) if path.exists() else default


def build_reviewer_package(output, prefix=None):
    from openpyxl import Workbook
    from openpyxl.styles import Font,PatternFill,Alignment
    from openpyxl.utils import get_column_letter
    root=Path(output);source=root/'genbank_submission';destination=root/'final_annotation_package'
    destination.mkdir(exist_ok=True)
    record=SeqIO.read(source/'genome.gbk','genbank');prefix=prefix or record.id
    prefix=''.join(c if c.isalnum() or c in '._-' else '_' for c in prefix)
    files={source/'genome.gbk':prefix+'.gbk',root/'genes.gff3':prefix+'.gff3',source/'features.tbl':prefix+'.tbl',source/'genome.fsa':prefix+'.fsa',root/'proteins.faa':prefix+'.faa',root/'cds.fna':prefix+'.ffn',root/'phagemine_qc_report.md':'phagemine_qc_report.md',root/'phagemine_qc_report.json':'phagemine_qc_report.json'}
    for src,name in files.items():shutil.copy2(src,destination/name)
    proteins={p['protein_id']:p for p in read(root/'evidence.json',[])}
    metadata=read(source/'submission_metadata.json',{})
    if metadata.get('locus_tag_prefix'):
        proteins.update({metadata['locus_tag_prefix'].rstrip('_')+'_'+p['protein_id'].rsplit('_',1)[-1]:p for p in list(proteins.values())})
    candidate_scores={(r['start'],r['end'],r['strand']):r for r in read(root/'gene_calls/start_refinement/start_candidates.json',[])+read(root/'gene_calls/start_refinement/short_orf_candidates.json',[])}
    cds=[f for f in record.features if f.type=='CDS'];qc=read(root/'phagemine_qc_report.json',{})
    wb=Workbook();wb.remove(wb.active)
    def sheet(name,columns,rows):
        ws=wb.create_sheet(name);ws.append(columns)
        for row in rows:
            values=[row.get(c) for c in columns] if isinstance(row,dict) else row
            ws.append([json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for v in values])
        for cells in ws:
            for cell in cells:
                if isinstance(cell.value,str):cell.data_type='s'  # annotation text is never an Excel formula
                cell.alignment=Alignment(vertical='top',wrap_text=True)
        for cell in ws[1]:cell.font=Font(bold=True,color='FFFFFF');cell.fill=PatternFill('solid',fgColor='17365D')
        ws.freeze_panes='B2';ws.auto_filter.ref=ws.dimensions;ws.row_dimensions[1].height=42
        for i,c in enumerate(columns,1):ws.column_dimensions[get_column_letter(i)].width=48 if any(w in c for w in ('Sequence','Translation','Flags','Reason','Hits')) else max(14,min(28,len(c)+2))
        return ws
    rows=[]
    for f in cds:
        q=f.qualifiers;tag=q['locus_tag'][0];pid=q.get('protein_id',[tag])[0].split('|')[-1]
        p=proteins.get(tag,proteins.get(pid,{}));params=p.get('gene_call_parameters',{});confidence=params.get('structural_confidence',{})
        strand=f.location.strand;start=int(f.location.start)+1;end=int(f.location.end)
        p5=isinstance(f.location.start,BeforePosition) if strand==1 else isinstance(f.location.end,AfterPosition)
        p3=isinstance(f.location.end,AfterPosition) if strand==1 else isinstance(f.location.start,BeforePosition)
        nt=str(f.extract(record.seq));cs=int(q.get('codon_start',['1'])[0]);coding=nt[cs-1:];protein=q.get('translation',[''])[0]
        seq=str(record.seq if strand==1 else record.seq.reverse_complement());left=start-1 if strand==1 else len(record)-end;rbs=rbs_score(seq,left)
        event=params.get('start_refinement',{});scores=candidate_scores.get((start,end,'+' if strand==1 else '-'),event.get('score_components',{}))
        errors=[]
        if cs not in (1,2,3) or (not p5 and cs!=1): errors.append('INVALID_CODON_START')
        if (p5 and (start!=1 if strand==1 else end!=len(record))) or (p3 and (end!=len(record) if strand==1 else start!=1)): errors.append('PARTIAL_NOT_AT_EDGE')
        usable=coding[:len(coding)//3*3];translated=str(Seq(usable).translate(table=11));translated=translated[:-1] if translated.endswith('*') else translated
        if not p5 and coding[:3] in {'ATG','GTG','TTG','CTG','ATT','ATC','ATA'} and translated:translated='M'+translated[1:]
        if translated!=protein:errors.append('TRANSLATION_MISMATCH')
        if not p5 and coding[:3] not in {'ATG','GTG','TTG','CTG','ATT','ATC','ATA'}:errors.append('INVALID_START')
        if not p3 and coding[-3:] not in {'TAA','TAG','TGA'}:errors.append('MISSING_STOP')
        if '*' in translated:errors.append('INTERNAL_STOP')
        if len(coding)%3 and not p3:errors.append('FRAME_ERROR')
        row=dict(zip(COLUMNS,[None]*len(COLUMNS)))
        row.update(Locus_Tag=tag,Start=start,End=end,Strand='+' if strand==1 else '-',Length_nt=len(nt),Length_aa=len(protein),Partial_5_Prime=p5,Partial_3_Prime=p3,Codon_Start=cs,Start_Codon=coding[:3],Stop_Codon=coding[-3:],Rare_Start_Flag=coding[:3] in {'CTG','ATT','ATC','ATA'},Upstream_20bp_DNA=seq[max(0,left-20):left],RBS_Motif=rbs['motif'],RBS_Spacer_nt=rbs['spacer_nt'],RBS_Score=scores.get('rbs_score',rbs['score']),Start_Total_Score=scores.get('score'),Start_Selection_Reason=params.get('overlap_rescue',{}).get('code') or event.get('code') or params.get('start_selection',{}).get('rule') or 'ORIGINAL_CALLER_START',Caller_Support=confidence.get('caller_support',p.get('gene_call_source')),Confidence_Tier=confidence.get('confidence_tier','unassessed'),Adjacent_Gene_Overlap_nt=max([max(0,min(end,int(g.location.end))-max(start,int(g.location.start)+1)+1) for g in cds if g is not f] or [0]),Product_Name=q.get('product',[''])[0],Codon_Verifier_Status='PASS' if not errors else '; '.join(errors),Review_Flags=params.get('start_selection_review',''),Protein_Translation=protein,Nucleotide_Sequence=nt)
        evidence=[e for e in p.get('evidence',[]) if e.get('status','').upper()=='REAL' and e.get('supports')]
        for names,hitcol,evalcol in [({'phrogs'},'Top_PHROGs_Hit','PHROGs_Evalue'),({'vogdb'},'Top_VOGDB_Hit','VOGDB_Evalue'),({'swissprot','swiss-prot'},'Top_SwissProt_Hit','SwissProt_Evalue')]:
            hits=[]
            for e in evidence:
                if e.get('source','').lower() not in names:continue
                m=e.get('metrics',{});value=m.get('evalue',m.get('e_value',m.get('full_sequence_e_value')))
                if value is not None:hits.append((float(value),e.get('identifier')))
            if hits:row[evalcol],row[hitcol]=min(hits)
        row['Domain_Hits']=[{'source':e.get('source'),'id':e.get('identifier'),'metrics':e.get('metrics')} for e in evidence if e.get('source','').lower() in {'pfam','cdd','interpro'}]
        rows.append(row)
    sheet('CDS_Verification',COLUMNS,rows)
    rna=read(root/'rna_features.json',{}).get('features',[])
    sheet('tRNA_tmRNA',['feature_id','type','start','end','strand','product','anticodon','source'],rna)
    queue=[r for r in read(root/'gene_calls/start_refinement/short_orf_candidates.json',[]) if not r.get('accepted')]
    overlaps={(r['start'],r['end'],r['strand']):r.get('families',[]) for r in read(root/'gene_calls/start_refinement/overlap_rescue.json',{}).get('candidates',[])}
    for row in queue:
        row['homology']={'nterminal':row.get('homology'),'family_hits':overlaps.get((row['start'],row['end'],row['strand']),[])}
    for alert in read(root/'ribosomal_slippage_review.json',[]):
        models=[proteins[i] for i in alert['cds_ids'] if i in proteins]
        if models: alert.update(start=min(p['start'] for p in models),end=max(p['end'] for p in models),strand=models[0]['strand'])
        queue.append(alert)
    queue+=[{**r,'review_reason':'STRICT_PRECISION_EXCLUDED'} for r in read(root/'cds_confidence.json',[]) if r.get('excluded')]
    sheet('Review_Queue_Candidates',['protein_id','start','end','strand','rbs_score','homology','review_reason','code','cds_ids','target','status'],queue)
    summary=[['Genome',record.id],['Length_nt',len(record)],['GC_percent',100*sum(b in 'GC' for b in str(record.seq))/len(record)],['CDS_count',len(cds)],['RNA_count',len(rna)],['Start_codons',dict(collections.Counter(r['Start_Codon'] for r in rows))]]
    summary += [[f"Guide_{c['question']}: {c['criterion']}",{'status':c['status'],'reason':c['reason']}] for c in qc.get('criteria',[])]
    sheet('Genome_Summary_and_QC',['Metric','Value'],summary)
    target=destination/(prefix+'_manual_verification.xlsx');wb.save(target)
    return destination
