import importlib.util
from pathlib import Path
import sys

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.SeqRecord import SeqRecord


def runner():
    directory = Path(__file__).resolve().parents[1] / 'evaluation/v1.3_validation/scripts'
    spec = importlib.util.spec_from_file_location('fresh_runner', directory / 'run_fresh.py')
    module = importlib.util.module_from_spec(spec)
    sys.path.insert(0, str(directory))
    try: spec.loader.exec_module(module)
    finally: sys.path.remove(str(directory))
    return module


def test_fresh_inputs_do_not_leak_reference_annotations(tmp_path):
    source = tmp_path / 'reference'; source.mkdir()
    output = tmp_path / 'output'; output.mkdir()
    record = SeqRecord(Seq('ATG' * 30), id='SYNTHETIC.1', name='SYNTHETIC',
                       annotations={'molecule_type': 'DNA'})
    record.features = [SeqFeature(FeatureLocation(0, 30, strand=1), type='CDS',
                                  qualifiers={'product': ['reference answer marker'], 'gene': ['answerA']}),
                       SeqFeature(CompoundLocation([FeatureLocation(35, 40, strand=1),
                                                     FeatureLocation(45, 55, strand=1)]), type='CDS')]
    SeqIO.write(record, source / 'fixture.gb', 'genbank')
    genomes, references, excluded = runner().prepare(source, output)
    fasta = Path(genomes[0]['fasta']).read_text()
    assert 'answer' not in fasta
    assert fasta == '>SYNTHETIC.1\n' + str(record.seq) + '\n'
    assert references[0]['gene'] == 'answerA'
    assert references[0]['start'] == 1 and references[0]['end'] == 30
    assert len(references) == 1 and excluded == 1
