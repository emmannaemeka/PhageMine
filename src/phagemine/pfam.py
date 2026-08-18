"""Optional local Pfam/HMMER domain evidence adapter.

No Pfam database is downloaded or bundled. A user must supply a local HMM file.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from .evidence import EvidenceAdapter, EvidenceAdapterResult
from .models import Evidence, EvidenceLevel, Protein


class PfamHMMAdapter(EvidenceAdapter):
    name = "PfamHMMAdapter"

    def __init__(self, pfam_path: str | Path | None = None, hmmscan: str | None = None, evalue_threshold: float | None = None, coverage_threshold: float | None = None, trusted_cutoff: bool = False, database_version: str | None = None, evidence_pack_id: str | None = None, threshold_mode: str | None = None, threads: int = 1):
        self.pfam_path = Path(pfam_path) if pfam_path else None
        self.hmmscan = hmmscan or shutil.which("hmmscan")
        self.evalue_threshold = evalue_threshold
        self.coverage_threshold = coverage_threshold
        self.threshold_mode = (threshold_mode or ("GA" if trusted_cutoff else ("MANUAL" if evalue_threshold is not None or coverage_threshold is not None else "NONE"))).upper()
        if self.threshold_mode not in {"GA", "MANUAL", "NONE"}:
            raise ValueError("Pfam threshold mode must be GA, MANUAL, or NONE")
        self.trusted_cutoff = self.threshold_mode == "GA"
        self.database_version = database_version or "unknown"
        self.evidence_pack_id = evidence_pack_id
        self.threads = max(1, int(threads))

    def available(self) -> bool:
        executable_available = bool(self.hmmscan and (Path(self.hmmscan).exists() or shutil.which(self.hmmscan)))
        return bool(executable_available and self.pfam_path and self.pfam_path.exists())

    def version(self) -> str:
        if not self.hmmscan or not (Path(self.hmmscan).exists() or shutil.which(self.hmmscan)):
            return "unavailable"
        result = subprocess.run([self.hmmscan, "-h"], capture_output=True, text=True, check=False)
        match = re.search(r"HMMER\s+([0-9][^\s]*)", (result.stdout or result.stderr))
        return match.group(1) if match else "unknown"

    def provenance(self) -> dict[str, Any]:
        return {"adapter": self.name, "adapter_version": "0.1.0", "hmmscan": self.hmmscan, "hmmer_version": self.version(), "pfam_path": str(self.pfam_path) if self.pfam_path else None, "evidence_pack_id": self.evidence_pack_id, "pfam_version": self.database_version, "threshold_mode": self.threshold_mode, "trusted_cutoff": self.trusted_cutoff, "thresholds": {"mode": self.threshold_mode, "evalue": self.evalue_threshold, "coverage": self.coverage_threshold, "trusted_cutoff": self.trusted_cutoff}, "status": "REAL" if self.available() else "UNAVAILABLE"}

    def analyze(self, proteins: list[Protein]) -> EvidenceAdapterResult:
        provenance = self.provenance()
        if not self.available():
            return EvidenceAdapterResult(self.name, "UNAVAILABLE", provenance=provenance, message="Pfam/HMMER unavailable; no domain evidence was fabricated.")
        with tempfile.TemporaryDirectory(prefix="phagemine-pfam-") as temp:
            fasta = Path(temp) / "proteins.faa"
            domtblout = Path(temp) / "hmmscan.domtblout"
            fasta.write_text("".join(f">{p.protein_id}\n{p.sequence}\n" for p in proteins))
            command = self._search_command(domtblout, fasta)
            result = subprocess.run(command, capture_output=True, text=True, check=False)
            provenance["search_command"] = command
            provenance["search_parameters"] = {"domtblout": True}
            if result.returncode:
                return EvidenceAdapterResult(self.name, "UNAVAILABLE", provenance=provenance, message=f"hmmscan failed with exit code {result.returncode}: {(result.stderr or result.stdout).strip()}")
            evidence = self.parse_domtblout(domtblout.read_text(), proteins, provenance)
        return EvidenceAdapterResult(self.name, "REAL", evidence=evidence, provenance=provenance)

    def _search_command(self, domtblout: Path, fasta: Path) -> list[str]:
        command = [str(self.hmmscan), "--cpu", str(self.threads), "--domtblout", str(domtblout), str(self.pfam_path), str(fasta)]
        if self.threshold_mode == "GA":
            command.insert(1, "--cut_ga")
        return command

    def _passes(self, evalue: float, coverage: float | None) -> tuple[bool, str]:
        if self.threshold_mode == "GA":
            return True, "STRONG"
        if self.threshold_mode == "NONE":
            return True, "WEAK"
        if self.evalue_threshold is not None and evalue > self.evalue_threshold:
            return False, "REJECTED"
        if self.coverage_threshold is not None and (coverage is None or coverage < self.coverage_threshold):
            return False, "REJECTED"
        if self.trusted_cutoff:
            return True, "STRONG"
        if self.evalue_threshold is not None and (self.coverage_threshold is None or coverage is not None):
            return True, "STRONG" if self.coverage_threshold is not None else "WEAK"
        return True, "WEAK"

    def parse_domtblout(self, text: str, proteins: list[Protein], provenance: dict[str, Any] | None = None) -> list[Evidence]:
        lengths = {protein.protein_id: protein.length for protein in proteins}
        records: list[Evidence] = []
        for line in text.splitlines():
            if not line or line.startswith("#"):
                continue
            fields = re.split(r"\s+", line.strip(), maxsplit=22)
            if len(fields) < 22:
                continue
            try:
                target, accession = fields[0], fields[1]
                protein_id = fields[3]
                full_evalue = float(fields[6])
                # HMMER 3.4 domtblout positions are retained explicitly. The
                # user-facing labels follow the C6 scan convention: independent
                # domain E-value is column 12 in the reported row, conditional is 13.
                independent_domain_evalue = float(fields[11])
                conditional_domain_evalue = float(fields[12])
                domain_evalue = independent_domain_evalue
                bit_score = float(fields[13])
                hmm_from, hmm_to = int(fields[15]), int(fields[16])
                ali_from, ali_to = int(fields[17]), int(fields[18])
                query_length = int(fields[5])
            except (ValueError, IndexError):
                continue
            if protein_id not in lengths:
                continue
            coverage = (ali_to - ali_from + 1) / query_length if query_length else None
            hmm_coverage = (hmm_to - hmm_from + 1) / int(fields[2]) if int(fields[2]) else None
            accepted, strength = self._passes(domain_evalue, coverage)
            description = fields[22] if len(fields) > 22 else target
            metrics = {"target_name": target, "pfam_accession": accession, "target_hmm_length": int(fields[2]), "query_protein_id": protein_id, "query_length": query_length, "full_sequence_e_value": full_evalue, "full_sequence_score": float(fields[7]), "conditional_domain_e_value": conditional_domain_evalue, "independent_domain_e_value": independent_domain_evalue, "domain_bit_score": bit_score, "hmm_coordinates": {"start": hmm_from, "end": hmm_to}, "query_coordinates": {"start": ali_from, "end": ali_to}, "envelope_coordinates": {"start": int(fields[19]), "end": int(fields[20])}, "accuracy": float(fields[21]), "query_coverage": coverage, "domain_coverage": hmm_coverage, "raw_hmmscan_row": line}
            level = EvidenceLevel.COMPUTATIONAL if strength == "STRONG" else EvidenceLevel.WEAK
            records.append(Evidence("domain", f"Pfam domain evidence for {target} ({accession}); this is evidence, not an automatic functional assignment.", level, "Pfam", self.database_version, status="REAL", supports=accepted, metrics=metrics, identifier=accession, threshold={"evalue": self.evalue_threshold, "coverage": self.coverage_threshold, "trusted_cutoff": self.trusted_cutoff}, coordinates={"start": ali_from, "end": ali_to}, provenance={**(provenance or {}), "protein_id": protein_id}, evidence_strength=strength.upper(), family_name=target, description=description))
        return records
