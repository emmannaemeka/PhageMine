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

    def __init__(self, pfam_path: str | Path | None = None, hmmscan: str | None = None, evalue_threshold: float | None = None, coverage_threshold: float | None = None, trusted_cutoff: bool = False, database_version: str | None = None, evidence_pack_id: str | None = None):
        self.pfam_path = Path(pfam_path) if pfam_path else None
        self.hmmscan = hmmscan or shutil.which("hmmscan")
        self.evalue_threshold = evalue_threshold
        self.coverage_threshold = coverage_threshold
        self.trusted_cutoff = trusted_cutoff
        self.database_version = database_version or "unknown"
        self.evidence_pack_id = evidence_pack_id

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
        return {"adapter": self.name, "adapter_version": "0.1.0", "hmmscan": self.hmmscan, "hmmer_version": self.version(), "pfam_path": str(self.pfam_path) if self.pfam_path else None, "evidence_pack_id": self.evidence_pack_id, "pfam_version": self.database_version, "thresholds": {"evalue": self.evalue_threshold, "coverage": self.coverage_threshold, "trusted_cutoff": self.trusted_cutoff}, "status": "REAL" if self.available() else "UNAVAILABLE"}

    def analyze(self, proteins: list[Protein]) -> EvidenceAdapterResult:
        provenance = self.provenance()
        if not self.available():
            return EvidenceAdapterResult(self.name, "UNAVAILABLE", provenance=provenance, message="Pfam/HMMER unavailable; no domain evidence was fabricated.")
        with tempfile.TemporaryDirectory(prefix="phagemine-pfam-") as temp:
            fasta = Path(temp) / "proteins.faa"
            domtblout = Path(temp) / "hmmscan.domtblout"
            fasta.write_text("".join(f">{p.protein_id}\n{p.sequence}\n" for p in proteins))
            command = [str(self.hmmscan), "--domtblout", str(domtblout), str(self.pfam_path), str(fasta)]
            if self.trusted_cutoff:
                command.insert(1, "--cut_ga")
            result = subprocess.run(command, capture_output=True, text=True, check=False)
            provenance["search_command"] = command
            provenance["search_parameters"] = {"domtblout": True}
            if result.returncode:
                return EvidenceAdapterResult(self.name, "UNAVAILABLE", provenance=provenance, message=f"hmmscan failed with exit code {result.returncode}: {(result.stderr or result.stdout).strip()}")
            evidence = self.parse_domtblout(domtblout.read_text(), proteins, provenance)
        return EvidenceAdapterResult(self.name, "REAL", evidence=evidence, provenance=provenance)

    def _passes(self, evalue: float, coverage: float | None) -> tuple[bool, str]:
        if self.evalue_threshold is not None and evalue > self.evalue_threshold:
            return False, "filtered_evalue"
        if self.coverage_threshold is not None and (coverage is None or coverage < self.coverage_threshold):
            return False, "filtered_coverage"
        if self.trusted_cutoff:
            return True, "strong"
        if self.evalue_threshold is not None and (self.coverage_threshold is None or coverage is not None):
            return True, "strong" if self.coverage_threshold is not None else "weak"
        return True, "weak"

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
                domain_evalue = float(fields[12])
                bit_score = float(fields[13])
                hmm_from, hmm_to = int(fields[15]), int(fields[16])
                ali_from, ali_to = int(fields[17]), int(fields[18])
                query_length = int(fields[5])
            except (ValueError, IndexError):
                continue
            if protein_id not in lengths:
                continue
            coverage = (ali_to - ali_from + 1) / query_length if query_length else None
            accepted, strength = self._passes(domain_evalue, coverage)
            if not accepted:
                continue
            description = fields[22] if len(fields) > 22 else target
            metrics = {"hmm_score": bit_score, "e_value": domain_evalue, "full_sequence_e_value": full_evalue, "independent_e_value": domain_evalue, "query_coverage": coverage, "hmm_coordinates": {"start": hmm_from, "end": hmm_to}, "domain_coordinates": {"start": ali_from, "end": ali_to}, "query_length": query_length}
            records.append(Evidence("domain", f"Pfam domain evidence for {target} ({accession}); this is evidence, not an automatic functional assignment.", EvidenceLevel.COMPUTATIONAL if strength == "strong" else EvidenceLevel.WEAK, "Pfam", self.database_version, status="REAL", supports=True, metrics=metrics, identifier=accession, threshold={"evalue": self.evalue_threshold, "coverage": self.coverage_threshold, "trusted_cutoff": self.trusted_cutoff}, coordinates={"start": ali_from, "end": ali_to}, provenance={**(provenance or {}), "protein_id": protein_id}, evidence_strength=strength, family_name=target, description=description))
        return records
