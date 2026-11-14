"""One-row-per-CDS view of the exact models supplied to GenBank export."""
import csv
from .genbank import _locus_tag, product_name
from .start_selection import STARTS_11, STOPS

COLUMNS = ['locus_tag', 'contig', 'start', 'end', 'strand', 'length_nt', 'length_aa',
           'start_codon', 'stop_codon', 'partial_5prime', 'partial_3prime', 'product',
           'review_flag', 'nucleotide_sequence', 'protein_sequence']


def write_orfs(path, proteins, overrides=None, locus_tag_prefix=None):
    with open(path, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter='\t')
        writer.writeheader()
        for p in proteins:
            product = (overrides or {}).get(p.protein_id, {}).get('product') or product_name(p)
            flags = []
            if p.cds[:3] in {"ATT", "ATC", "ATA", "CTG"}: flags.append("rare start codon - review")
            if p.gene_call_parameters.get("start_selection_review"):
                flags.append(p.gene_call_parameters["start_selection_review"])
            if p.partial_5prime or p.cds[:3] not in STARTS_11: flags.append('unresolved_5prime')
            if p.partial_3prime or p.cds[-3:] not in STOPS: flags.append('unresolved_3prime')
            if '*' in p.sequence: flags.append('internal_stop_review_no_automatic_recoding')
            if product.lower() == 'hypothetical protein': flags.append('no_supported_product')
            writer.writerow(dict(locus_tag=_locus_tag(p, locus_tag_prefix), contig=p.genome_id,
                                 start=p.start, end=p.end, strand=p.strand, length_nt=len(p.cds),
                                 length_aa=len(p.sequence), start_codon=p.cds[:3], stop_codon=p.cds[-3:],
                                 partial_5prime=p.partial_5prime, partial_3prime=p.partial_3prime,
                                 product=product, review_flag=';'.join(flags),
                                 nucleotide_sequence=p.cds, protein_sequence=p.sequence))
