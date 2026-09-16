"""Regenerate compact benchmark figures from the distributed TSV tables."""
from pathlib import Path
import csv
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1]; tab=root/'tables'; fig=root/'figures'
rows=list(csv.DictReader(open(tab/'strict_relaxed_metrics.tsv'),delimiter='\t'))
tools=['PhageMine','Pharokka','Prokka']; x=range(3)
def agg(t,k): return [float(r[k]) for r in rows if r['tool']==t]
fig1,ax=plt.subplots(figsize=(6,4)); strict=[]; relaxed=[]
for t in tools:
 a=[r for r in rows if r['tool']==t]; strict.append(2*sum(int(r['strict_tp']) for r in a)/(sum(int(r['predicted_cds'])+int(r['reference_cds']) for r in a))); relaxed.append(sum(float(r['relaxed_f1']) for r in a)/len(a))
ax.bar([i-.18 for i in x],strict,.36,label='Strict exact-coordinate F1'); ax.bar([i+.18 for i in x],relaxed,.36,label='Relaxed overlap F1'); ax.set_xticks(list(x),tools); ax.set_ylim(0,1); ax.set_ylabel('F1'); ax.legend(frameon=False); fig1.tight_layout(); fig1.savefig(fig/'strict_vs_relaxed_f1_regenerated.png',dpi=300)
