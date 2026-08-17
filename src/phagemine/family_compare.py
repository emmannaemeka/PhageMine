"""Descriptive, offline analytics for a persisted PMF database."""
from __future__ import annotations
import csv, json
from collections import Counter, defaultdict
from pathlib import Path

def load_families(root):
    root=Path(root); p=root/'families.json' if root.is_dir() else root
    families=json.loads(p.read_text()); members=[]
    mp=root/'family_members.tsv'
    if mp.exists():
        with mp.open() as h: members=list(csv.DictReader(h,delimiter='\t'))
    return families,members

def _genomes(families,members):
    return sorted({m.get('source_genome_id','UNKNOWN') for m in members})

def compare_database(root, output, core_genomes=None):
    families,members=load_families(root); out=Path(output); out.mkdir(parents=True,exist_ok=True)
    byfam=defaultdict(list)
    for m in members: byfam[m['family_id']].append(m)
    genomes=_genomes(families,members); famids=sorted(byfam)
    presence={g:Counter() for g in genomes}
    for fid,ms in byfam.items():
        for g,n in Counter(m.get('source_genome_id','UNKNOWN') for m in ms).items(): presence[g][fid]=n
    stats={'total_proteins':len(members),'total_families':len(famids),
           'singleton_families':sum(len(v)==1 for v in byfam.values()),
           'multi_member_families':sum(len(v)>1 for v in byfam.values()),
           'genome_count':len(genomes),'host_genera':sorted({m['host_genus'] for m in members if m.get('host_genus')}),
           'family_size_distribution':dict(sorted(Counter(len(v) for v in byfam.values()).items()))}
    pair=[]
    for i,a in enumerate(genomes):
        for b in genomes[i+1:]:
            A={f for f in famids if presence[a][f]}; B={f for f in famids if presence[b][f]}; u=A|B
            pair.append({'genome_a':a,'genome_b':b,'proteins_a':sum(presence[a].values()),'proteins_b':sum(presence[b].values()),'shared_pmf_count':len(A&B),'private_a_count':len(A-B),'private_b_count':len(B-A),'jaccard_similarity':len(A&B)/len(u) if u else None})
    selected=core_genomes or genomes; core=[f for f in famids if all(presence.get(g,Counter())[f] for g in selected)]
    with (out/'pmf_presence_absence.tsv').open('w',newline='') as h:
        w=csv.writer(h,delimiter='\t'); w.writerow(['family_id']+genomes)
        for f in famids:w.writerow([f]+[presence[g][f] for g in genomes])
    with (out/'pmf_pairwise_genomes.tsv').open('w',newline='') as h:
        w=csv.DictWriter(h,fieldnames=pair[0].keys() if pair else ['genome_a','genome_b'],delimiter='\t');w.writeheader();w.writerows(pair)
    with (out/'pmf_core_families.tsv').open('w') as h: h.write('family_id\n'+'\n'.join(core)+('\n' if core else ''))
    hosts={}
    for f,ms in byfam.items():
        hs={m.get('host_genus') for m in ms if m.get('host_genus')};
        if not hs: state='HOST_METADATA_INCOMPLETE'
        elif hs=={'Pseudomonas'}: state='OBSERVED_IN_PSEUDOMONAS_COHORT_ONLY'
        elif hs=={'Salmonella'}: state='OBSERVED_IN_SALMONELLA_COHORT_ONLY'
        elif hs & {'Pseudomonas','Salmonella'}: state='CROSS_HOST_COHORT'
        else: state='HOST_METADATA_INCOMPLETE'
        hosts[f]=state
    with (out/'pmf_host_group_summary.tsv').open('w',newline='') as h:
        w=csv.writer(h,delimiter='\t');w.writerow(['family_id','host_group_observation']);w.writerows(hosts.items())
    summary={'statistics':stats,'genomes':genomes,'core_genomes':selected,'core_family_count':len(core),'pairwise':pair,'host_group_counts':dict(Counter(hosts.values()))}
    (out/'pmf_comparison_summary.json').write_text(json.dumps(summary,indent=2,sort_keys=True))
    with (out/'pmf_comparison_summary.tsv').open('w') as h:
        for k,v in stats.items(): h.write(f'{k}\t{json.dumps(v,sort_keys=True)}\n')
    return summary

def sensitivity_rows(builds, core_genomes=None):
    rows=[]
    for identity,summary in sorted(builds.items()):
        stats=summary['statistics']; rows.append({'minimum_identity':identity,'minimum_coverage':summary.get('minimum_coverage'), 'coverage_mode':summary.get('coverage_mode',0),'clustering_mode':summary.get('clustering_mode',0),'backend':summary.get('backend'),'total_families':stats['total_families'],'singleton_families':stats['singleton_families'],'multi_member_families':stats['multi_member_families'],'largest_family_size':max(map(int,stats['family_size_distribution']) if stats['family_size_distribution'] else [0]),'core_all_genomes':summary.get('core_family_count',0),'cross_host_cohort_families':summary.get('host_group_counts',{}).get('CROSS_HOST_COHORT',0)})
    return rows
