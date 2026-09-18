"""Plot adjudication status, explicitly retaining unresolved cases."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]

def main():
    with (ROOT/'figure_source.tsv').open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    plt.rcParams.update({'svg.hashsalt':'phagemine-functional-v1','pdf.fonttype':42})
    fig,axes=plt.subplots(1,2,figsize=(10,4.8),sharey=True)
    groups=[('Exact / explicit synonym',{'EXACT_PRODUCT_AGREEMENT','EQUIVALENT_FUNCTION'},'#0072B2'),
            ('Compatible (reviewed)',{'COMPATIBLE_BUT_BROADER','COMPATIBLE_BUT_MORE_SPECIFIC'},'#009E73'),
            ('Unsupported / disagreement',{'UNSUPPORTED_SPECIFICITY','GENUINE_FUNCTIONAL_DISAGREEMENT'},'#D55E00'),
            ('Abstention',{'ABSTENTION'},'#999999'),
            ('Unresolved — review required',{'UNRESOLVED_REVIEW_REQUIRED'},'#E69F00')]
    for ax,mode in zip(axes,['exact','relaxed']):
        bottom=[0,0]
        for label,categories,color in groups:
            values=[sum(int(r['count']) for r in rows if r['matching_mode']==mode and r['tool']==t and r['category'] in categories) for t in ['PhageMine','Pharokka']]
            bars=ax.bar(['PhageMine','Pharokka'],values,bottom=bottom,label=label,color=color)
            for bar,value,base in zip(bars,values,bottom):
                if value: ax.text(bar.get_x()+bar.get_width()/2,base+value/2,str(value),ha='center',va='center',fontsize=9)
            bottom=[a+b for a,b in zip(bottom,values)]
        ax.set_title('Exact loci' if mode=='exact' else 'Relaxed, unambiguous loci')
        ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel('Matched, evaluable named-reference loci')
    fig.suptitle('Reference concordance: manual adjudication incomplete')
    handles,labels=axes[0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',ncol=2,frameon=False,fontsize=9)
    fig.tight_layout(rect=(0,.19,1,.95))
    (ROOT/'figures').mkdir(exist_ok=True)
    for ext in ['png','svg','pdf']:
        metadata={'Date':None} if ext=='svg' else {'CreationDate':None,'ModDate':None} if ext=='pdf' else {}
        fig.savefig(ROOT/'figures'/f'adjudication_status.{ext}',dpi=300,metadata=metadata)
    svg = ROOT / 'figures/adjudication_status.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    plt.close(fig)

if __name__=='__main__': main()
