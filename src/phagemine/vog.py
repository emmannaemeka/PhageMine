"""Optional local VOGDB viral-orthology evidence adapter."""
from __future__ import annotations

import gzip
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from .evidence import EvidenceAdapter, EvidenceAdapterResult
from .models import Evidence, EvidenceLevel, Protein


class VOGHMMAdapter(EvidenceAdapter):
    name = "VOGHMMAdapter"

    def __init__(self, vog_path: str | Path | None = None, annotations_path: str | Path | None = None,
                 hmmscan: str | None = None, evalue_threshold: float | None = 1e-5,
                 coverage_threshold: float | None = 0.5, database_version: str | None = None):
        self.vog_path = Path(vog_path).expanduser() if vog_path else None
        self.annotations_path = Path(annotations_path).expanduser() if annotations_path else None
        self.hmmscan = hmmscan or shutil.which("hmmscan")
        self.evalue_threshold = evalue_threshold
        self.coverage_threshold = coverage_threshold
        self.database_version = database_version or "unknown"

    def available(self) -> bool:
        executable = bool(self.hmmscan and (Path(self.hmmscan).exists() or shutil.which(self.hmmscan)))
        return bool(executable and self.vog_path and self.vog_path.is_file())

    def version(self) -> str:
        if not self.hmmscan or not (Path(self.hmmscan).exists() or shutil.which(self.hmmscan)):
            return "unavailable"
        result = subprocess.run([self.hmmscan, "-h"], capture_output=True, text=True, check=False)
        match = re.search(r"HMMER\s+([0-9][^\s]*)", (result.stdout or result.stderr))
        return match.group(1) if match else "unknown"

    def provenance(self) -> dict[str, Any]:
        return {"adapter": self.name, "adapter_version": "0.1.0", "hmmscan": self.hmmscan,
                "hmmer_version": self.version(), "vog_path": str(self.vog_path) if self.vog_path else None,
                "annotations_path": str(self.annotations_path) if self.annotations_path else None,
                "vogdb_version": self.database_version, "threshold_mode": "MANUAL",
                "thresholds": {"evalue": self.evalue_threshold, "coverage": self.coverage_threshold},
                "trusted_cutoff": False, "status": "REAL" if self.available() else "UNAVAILABLE"}

    def analyze(self, proteins: list[Protein]) -> EvidenceAdapterResult:
        provenance = self.provenance()
        if not self.available():
            return EvidenceAdapterResult(self.name, "UNAVAILABLE", provenance=provenance,
                                         message="VOGDB/HMMER unavailable or not a prepared combined HMM; no VOG evidence was fabricated.")
        with tempfile.TemporaryDirectory(prefix="phagemine-vog-") as temp:
            fasta = Path(temp) / "proteins.faa"
            domtblout = Path(temp) / "hmmscan.domtblout"
            fasta.write_text("".join(f">{p.protein_id}\n{p.sequence}\n" for p in proteins))
            command = [str(self.hmmscan), "--domtblout", str(domtblout), str(self.vog_path), str(fasta)]
            result = subprocess.run(command, capture_output=True, text=True, check=False)
            provenance["search_command"] = command
            provenance["search_parameters"] = {"domtblout": True, "threshold_mode": "MANUAL"}
            if result.returncode:
                return EvidenceAdapterResult(self.name, "UNAVAILABLE", provenance=provenance,
                    message=f"hmmscan failed with exit code {result.returncode}: {(result.stderr or result.stdout).strip()}")
            evidence = self.parse_domtblout(domtblout.read_text(), proteins, provenance)
        return EvidenceAdapterResult(self.name, "REAL", evidence=evidence, provenance=provenance)

    def _passes(self, evalue: float, coverage: float | None) -> tuple[bool, str]:
        if self.evalue_threshold is not None and evalue > self.evalue_threshold:
            return False, "REJECTED"
        if self.coverage_threshold is not None and (coverage is None or coverage < self.coverage_threshold):
            return False, "REJECTED"
        if self.evalue_threshold is None and self.coverage_threshold is None:
            return True, "WEAK"
        return True, "STRONG"

    def _annotations(self) -> dict[str, dict[str, str | None]]:
        if not self.annotations_path or not self.annotations_path.exists():
            return {}
        opener = gzip.open if self.annotations_path.suffix == ".gz" else open
        with opener(self.annotations_path, "rt", encoding="utf-8") as handle:
            header = [item.lstrip("#") for item in handle.readline().rstrip("\n").split("\t")]
            result = {}
            for line in handle:
                values = line.rstrip("\n").split("\t")
                if len(values) < len(header):
                    continue
                row = dict(zip(header, values))
                vog_id = row.get("GroupName") or row.get("VOG") or row.get("VOG_ID")
                if vog_id:
                    result[vog_id] = {"functional_category": row.get("FunctionalCategory") or None,
                                      "description": row.get("ConsensusFunctionalDescription") or None}
            return result

    def parse_domtblout(self, text: str, proteins: list[Protein], provenance: dict[str, Any] | None = None) -> list[Evidence]:
        lengths = {p.protein_id: p.length for p in proteins}
        annotations = self._annotations()
        records: list[Evidence] = []
        for line in text.splitlines():
            if not line or line.startswith("#"):
                continue
            fields = re.split(r"\s+", line.strip(), maxsplit=22)
            if len(fields) < 22:
                continue
            try:
                vog_id, protein_id = fields[0], fields[3]
                target_length, query_length = int(fields[2]), int(fields[5])
                full_evalue, full_score = float(fields[6]), float(fields[7])
                independent = float(fields[11]); conditional = float(fields[12]); bit_score = float(fields[13])
                hmm_from, hmm_to = int(fields[15]), int(fields[16])
                query_from, query_to = int(fields[17]), int(fields[18])
                env_from, env_to = int(fields[19]), int(fields[20]); accuracy = float(fields[21])
            except (ValueError, IndexError):
                continue
            if protein_id not in lengths:
                continue
            coverage = (query_to - query_from + 1) / query_length if query_length else None
            accepted, strength = self._passes(independent, coverage)
            ann = annotations.get(vog_id, {})
            hmmscan_description = fields[22] if len(fields) > 22 else None
            description = ann.get("description") or hmmscan_description
            metrics = {"vog_id": vog_id, "query_protein_id": protein_id, "query_length": query_length,
                       "full_sequence_e_value": full_evalue, "full_sequence_score": full_score,
                       "conditional_domain_e_value": conditional, "independent_domain_e_value": independent,
                       "domain_bit_score": bit_score, "hmm_coordinates": {"start": hmm_from, "end": hmm_to},
                       "query_coordinates": {"start": query_from, "end": query_to},
                       "envelope_coordinates": {"start": env_from, "end": env_to}, "accuracy": accuracy,
                       "query_coverage": coverage, "domain_coverage": (hmm_to - hmm_from + 1) / target_length if target_length else None,
                       "raw_hmmscan_row": line, "hmmscan_description": hmmscan_description,
                       "consensus_functional_description": ann.get("description")}
            records.append(Evidence("viral_orthology", f"VOGDB orthologous-group evidence for {vog_id}; evidence only, not an automatic functional assignment.",
                EvidenceLevel.COMPUTATIONAL if strength == "STRONG" else EvidenceLevel.WEAK, "VOGDB", self.database_version,
                status="REAL", supports=(accepted and strength != "REJECTED"), metrics=metrics, identifier=vog_id,
                threshold={"evalue": self.evalue_threshold, "coverage": self.coverage_threshold, "threshold_mode": "MANUAL"},
                coordinates={"start": query_from, "end": query_to}, provenance={**(provenance or {}), "protein_id": protein_id, "threshold_mode": "MANUAL", "vogdb_version": self.database_version},
                evidence_strength=strength, family_name=vog_id, description=description))
            records[-1].metrics["functional_category"] = ann.get("functional_category")
            records[-1].metrics["functional_category_meaning"] = None
        return records
