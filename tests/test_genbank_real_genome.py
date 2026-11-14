"""Export invariants using real sequence and independently saved PHANOTATE calls."""
import csv
import json
from pathlib import Path
from Bio import SeqIO
from phagemine.models import Protein, SubmissionMetadata, EvidenceLevel
from phagemine.genbank import write_package
from phagemine.genbank_flatfile import write_flatfile
from phagemine.reporting import write_final_gene_models

FIXTURE = Path(__file__).parent / 'fixtures/genbank_real'

def test_real_genome_export_invariants(tmp_path):
    genome = SeqIO.read(FIXTURE / 'PZ583544.1.fasta', 'fasta')
    proteins = []
    with (FIXTURE / 'PZ583544.1.models.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            start, end = int(row['start']), int(row['end'])
            cds = genome.seq[start-1:end]
            if row['strand'] == '-':
                cds = cds.reverse_complement()
            proteins.append(Protein(genome.id, row['protein_id'], start, end,
                                    row['strand'], str(cds), row['sequence'], 'PHANOTATE',
                                    annotation='supported computational regression label',
                                    annotation_level=EvidenceLevel.COMPUTATIONAL))
    write_final_gene_models(tmp_path, proteins)
    write_package(tmp_path, genome.id, str(genome.seq), proteins, {'input_sha256': 'fixture'},
                  SubmissionMetadata(topology='linear'), table2asn_executable='not-installed')
    record = SeqIO.read(tmp_path / 'genbank_submission/genome.gbk', 'genbank')
    with (tmp_path / 'final_gene_models.tsv').open() as handle:
        models = list(csv.DictReader(handle, delimiter='\t'))
    genes = [f for f in record.features if f.type == 'gene']
    cds_features = [f for f in record.features if f.type == 'CDS']
    assert len(genes) == len(cds_features) == len(models) == 72
    assert len(record.name) <= 16
    assert record.annotations['topology'] == 'linear'
    provenance = json.loads((tmp_path / 'genbank_submission/cds_provenance.json').read_text())
    validation = json.loads((tmp_path / 'genbank_submission/validation.json').read_text())
    assert validation['phagemine_pre_submission_validation']['errors'] == []
    assert provenance['PM_000001']['partial_5prime']
    assert provenance['PM_000072']['partial_3prime']
    for gene, cds, model in zip(genes, cds_features, models):
        exported = provenance[model['final_protein_id']]
        if model['final_protein_id'] == 'PM_000072':
            assert int(cds.location.end) == len(genome) == int(model['end']) + 1
        else:
            assert int(cds.location.end) == int(model['end'])
        assert gene.location == cds.location
        assert int(cds.location.start) + 1 == int(model['start'])
        assert cds.location.strand == (1 if model['strand'] == '+' else -1)
        coding_sequence = cds.extract(record.seq)
        observed = str(coding_sequence[:len(coding_sequence) // 3 * 3].translate(table=11)).rstrip('*')
        assert observed == model['translation']
        if not exported['partial_5prime']:
            observed = 'M' + observed[1:]
        assert cds.qualifiers['translation'][0] == observed
        assert cds.qualifiers['product'] == ['supported computational regression label']
    package = tmp_path / 'genbank_submission'
    sequence = (package / 'genome.fsa').read_text()
    long_id = 'a_very_long_real_genome_identifier'
    (package / 'genome.fsa').write_text(sequence.replace('>'+genome.id, '>'+long_id, 1))
    write_flatfile(package)
    renamed = SeqIO.read(package / 'genome.gbk', 'genbank')
    assert len(renamed.name) <= 16 and renamed.id == long_id
