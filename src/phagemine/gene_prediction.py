"""Gene-prediction adapters. Production analysis selects an external predictor explicitly."""
from __future__ import annotations

import re
import shutil
import subprocess
from abc import ABC, abstractmethod
from pathlib import Path

from .genome import predict_orfs, translate
from .models import Protein


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def _assign_ids(proteins: list[Protein]) -> list[Protein]:
    proteins.sort(key=lambda item: (item.start, item.end, item.strand))
    for index, protein in enumerate(proteins, 1):
        protein.protein_id = f"PM_{index:06d}"
        protein.locus_tag = protein.protein_id
        protein.start_codon = protein.cds[:3] if len(protein.cds) >= 3 else None
        protein.stop_codon = protein.cds[-3:] if len(protein.cds) >= 3 else None
    return proteins


class GenePredictor(ABC):
    name: str

    @abstractmethod
    def predict(self, genome_id: str, sequence: str, input_fasta: str | Path | None = None) -> list[Protein]: ...

    @abstractmethod
    def version(self) -> str: ...

    @abstractmethod
    def parameters(self) -> dict: ...


class DemoORFPredictor(GenePredictor):
    """Deterministic test-fixture caller; never a production default."""
    name = "simple_orf_demo"

    def __init__(self, min_orf_nt: int = 90):
        self.min_orf_nt = min_orf_nt

    def predict(self, genome_id: str, sequence: str, input_fasta: str | Path | None = None) -> list[Protein]:
        return _assign_ids(predict_orfs(genome_id, sequence, self.min_orf_nt))

    def version(self) -> str:
        return "0.1.0"

    def parameters(self) -> dict:
        return {"min_orf_nt": self.min_orf_nt, "status": "test fixture only"}


class PHANOTATEPredictor(GenePredictor):
    """Adapter for the locally installed PHANOTATE command-line predictor."""
    name = "PHANOTATE"

    def __init__(self, executable: str | None = None, extra_args: list[str] | None = None):
        self.executable = executable or shutil.which("phanotate.py") or shutil.which("phanotate")
        self.extra_args = extra_args or []

    def available(self) -> bool:
        return bool(self.executable and Path(self.executable).exists())

    def version(self) -> str:
        if not self.available():
            return "unavailable"
        for flag in ("--version", "-v"):
            result = subprocess.run([str(self.executable), flag], capture_output=True, text=True, check=False)
            text = (result.stdout or result.stderr).strip()
            if text:
                return text.splitlines()[0]
        return "reported by executable at run time"

    def parameters(self) -> dict:
        return {"executable": self.executable, "extra_args": self.extra_args}

    def predict(self, genome_id: str, sequence: str, input_fasta: str | Path | None = None) -> list[Protein]:
        if not self.available():
            raise RuntimeError("PHANOTATE is required for production gene prediction but was not found. Install PHANOTATE locally and add phanotate.py (or phanotate) to PATH, or pass --phanotate /path/to/phanotate.py. The demo predictor is test-fixture only.")
        if input_fasta is None:
            raise ValueError("PHANOTATE requires the original FASTA path for reproducible external execution.")
        completed = subprocess.run([str(self.executable), *self.extra_args, str(input_fasta)], capture_output=True, text=True, check=False)
        if completed.returncode:
            raise RuntimeError(f"PHANOTATE failed (exit {completed.returncode}): {(completed.stderr or completed.stdout).strip()}")
        return self.parse_output(genome_id, sequence, completed.stdout)

    @staticmethod
    def parse_output(genome_id: str, sequence: str, output: str) -> list[Protein]:
        """Parse PHANOTATE's tab/whitespace coordinate table without assuming a fixed header."""
        proteins: list[Protein] = []
        for raw in output.splitlines():
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = re.split(r"\s+", line)
            if len(fields) < 3:
                continue
            try:
                raw_start, raw_end = int(fields[0]), int(fields[1])
            except ValueError:
                continue
            strand_token = fields[2]
            strand = "-" if strand_token in {"-", "-1"} or raw_start > raw_end else "+"
            start, end = min(raw_start, raw_end), max(raw_start, raw_end)
            if not (1 <= start <= end <= len(sequence)):
                raise ValueError(f"PHANOTATE emitted out-of-range coordinates: {raw_start}, {raw_end}")
            genomic_cds = sequence[start - 1:end]
            cds = reverse_complement(genomic_cds) if strand == "-" else genomic_cds
            protein = Protein(genome_id, "", start, end, strand, cds, translate(cds), "PHANOTATE")
            protein.gene_call_parameters = {"raw_start": raw_start, "raw_end": raw_end, "reported_strand": strand_token}
            proteins.append(protein)
        if not proteins:
            raise ValueError("PHANOTATE produced no parseable gene calls; inspect its output format and version.")
        return _assign_ids(proteins)


def create_predictor(name: str = "phanotate", executable: str | None = None) -> GenePredictor:
    if name.lower() == "phanotate":
        return PHANOTATEPredictor(executable)
    if name.lower() == "demo":
        return DemoORFPredictor()
    raise ValueError(f"Unknown gene predictor: {name}")
