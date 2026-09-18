"""Reconstruct published lexical yield directly from archived frozen annotations."""
from pathlib import Path
import re
from run_benchmark import REPO, digest, gff, read, write


def main():
    root=Path(__file__).resolve().parents[1]
    table=REPO/'docs/benchmark_v1.2/tables'
    rows=[]
    for source in read(table/'strict_relaxed_metrics.tsv'):
        tool,a=source['tool'],source['accession']
        path=root/'predictions'/tool/a/('annotation.tsv' if tool=='PhageMine' else f'{a}.gff')
        predictions=read(path) if tool=='PhageMine' else gff(path,a)
        total=len(predictions)
        assert total==int(source['predicted_cds'])
        named=sum(bool(p['product']) and not re.search('hypothetical|unknown function|uncharacterized|uncharacterised',p['product'],re.I) for p in predictions)
        rows.append(dict(genome=source['genome'],accession=a,tool=tool,predicted_cds=total,
            named_product_calls=named,hypothetical_or_uncharacterized=total-named,named_fraction=f'{named/total:.6f}',
            measurement_status='PER_GENOME_AVAILABLE',source_file=str(path.relative_to(REPO)),source_sha256=digest(path)))
    totals=[]
    for tool,expected,predicted in [('PhageMine',439,804),('Pharokka',398,804),('Prokka',328,675)]:
        subset=[r for r in rows if r['tool']==tool]
        n=sum(r['named_product_calls'] for r in subset); p=sum(r['predicted_cds'] for r in subset)
        assert (n,p)==(expected,predicted), (tool,n,p)
        totals.append(dict(tool=tool,predicted_cds=p,named_product_calls=n,hypothetical_or_uncharacterized=p-n,named_fraction=f'{n/p:.6f}'))
    write(table/'functional_yield_by_genome.tsv',rows)
    write(table/'functional_yield_totals.tsv',totals)

if __name__=='__main__': main()
