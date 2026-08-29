from pathlib import Path

from phagemine.gene_models import GeneModel
from phagemine.models import Protein
from phagemine.pipeline import _write_gene_call_provenance
from phagemine.genome_representation import GenomeRepresentation


def test_gene_model_retains_raw_identity_and_serializes_lengths():
    model = GeneModel(
        "PHANOTATE", "P1", 10, 39, "+", cds_sequence="ATG" * 10,
        protein_sequence="M" * 10, raw_start=10, raw_end=39,
        raw_strand="+", source_file="phanotate.raw.txt",
    )
    record = model.to_dict()
    assert model.raw_identifier == "P1"
    assert record["length_nt"] == 30
    assert record["length_aa"] == 10
    assert record["raw_start"] == 10
    assert record["coordinate_system"] == "1-based-inclusive"


def test_gene_call_provenance_writes_raw_stream_manifest_and_trace(tmp_path):
    class Predictor:
        name = "PHANOTATE"
        last_raw_output = "1 30 + genome\n"
        last_command = ["/bin/phanotate.py", "genome.fasta"]

        def version(self):
            return "1.6.7"

        def parameters(self):
            return {"executable": "/bin/phanotate.py", "extra_args": []}

    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">genome\n" + "ATG" * 10 + "\n")
    proteins = [Protein("genome", "PM_000001", 1, 30, "+", "ATG" * 10, "M" * 10, "PHANOTATE")]
    _write_gene_call_provenance(tmp_path / "out", fasta, GenomeRepresentation.original("genome", "ATG" * 10), Predictor(), proteins)
    root = tmp_path / "out" / "gene_calls"
    assert (root / "raw" / "phanotate.raw.txt").read_text() == "1 30 + genome\n"
    assert "raw_identifier" in (root / "raw" / "phanotate.tsv").read_text()
    assert (root / "gene_call_manifest.json").is_file()
    trace = (root / "final_gene_model_trace.tsv").read_text()
    assert "PM_000001" in trace
    assert "LEGACY_PHANOTATE_FINAL" in trace
