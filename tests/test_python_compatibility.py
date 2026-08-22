import ast
from pathlib import Path

from phagemine.genome_representation import GenomeRepresentation
from phagemine.models import Protein
from phagemine.reporting import write_outputs
from phagemine.sequencing_provenance import SequencingProvenance


def test_src_parses_as_python_3_10():
    root = Path(__file__).parents[1]
    for path in (root / "src").rglob("*.py"):
        ast.parse(path.read_text(), filename=str(path), feature_version=(3, 10))


def test_gff3_preserves_raw_coordinates(tmp_path):
    protein = Protein("genome", "PM_000001", 10, 20, "+", "ATG", "M", "fixture",
                      gene_call_parameters={"raw_start": 12, "raw_end": 19})
    representation = GenomeRepresentation.original("genome", "ATG")
    write_outputs(tmp_path, representation, SequencingProvenance(), [protein], [], {})
    gff = (tmp_path / "genes.gff3").read_text()
    assert "raw_start=12;raw_end=19" in gff
