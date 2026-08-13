"""Optional local PHROGs MMseqs2 profile-search evidence adapter."""
from __future__ import annotations

import csv
import gzip
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from .evidence import EvidenceAdapter, EvidenceAdapterResult
from .models import Evidence, EvidenceLevel, Protein


MMSEQS_FIELDS = (
    "query", "target", "fident", "alnlen", "qstart", "qend", "qlen",
    "tstart", "tend", "tlen", "evalue", "bits", "qcov", "tcov",
)
MMSEQS_FORMAT = ",".join(MMSEQS_FIELDS)


class PHROGSMMseqsAdapter(EvidenceAdapter):
    name = "PHROGSMMseqsAdapter"

    def __init__(self, database_path: str | Path | None = None,
                 annotations_path: str | Path | None = None,
                 mmseqs: str | None = None, database_version: str | None = None,
                 evalue_threshold: float | None = 1e-5,
                 coverage_threshold: float | None = 0.5,
                 score_threshold: float | None = None,
                 identity_threshold: float | None = None,
                 alignment_length_threshold: int | None = None):
        self.database_path = Path(database_path).expanduser() if database_path else None
        self.annotations_path = Path(annotations_path).expanduser() if annotations_path else None
        self.mmseqs = mmseqs or shutil.which("mmseqs")
        self.database_version = database_version or "unknown"
        self.evalue_threshold = evalue_threshold
        self.coverage_threshold = coverage_threshold
        self.score_threshold = score_threshold
        self.identity_threshold = identity_threshold
        self.alignment_length_threshold = alignment_length_threshold

    def available(self) -> bool:
        executable = bool(self.mmseqs and (Path(self.mmseqs).exists() or shutil.which(self.mmseqs)))
        # An MMseqs2 database is addressed by its prefix. Its dbtype companion
        # is the stable preparation marker; Pharokka also ships a _h database.
        prepared = bool(self.database_path and self.database_path.is_file() and
                        Path(f"{self.database_path}.dbtype").is_file())
        return bool(executable and prepared)

    def version(self) -> str:
        if not self.mmseqs or not (Path(self.mmseqs).exists() or shutil.which(self.mmseqs)):
            return "unavailable"
        result = subprocess.run([self.mmseqs, "version"], capture_output=True, text=True, check=False)
        output = (result.stdout or result.stderr).strip()
        match = re.search(r"(?:Version:\s*)?([0-9][0-9.\-A-Za-z]+)", output)
        return match.group(1) if match else (output.splitlines()[0] if output else "unknown")

    def provenance(self) -> dict[str, Any]:
        thresholds = self._thresholds()
        return {
            "adapter": self.name, "adapter_version": "0.2.0",
            "mmseqs": self.mmseqs, "mmseqs_version": self.version(),
            "phrogs_path": str(self.database_path) if self.database_path else None,
            "annotations_path": str(self.annotations_path) if self.annotations_path else None,
            "phrogs_version": self.database_version, "threshold_mode": "MANUAL",
            "output_format": MMSEQS_FORMAT,
            "thresholds": thresholds,
            "status": "REAL" if self.available() else "UNAVAILABLE",
        }

    def _thresholds(self) -> dict[str, Any]:
        return {"evalue": self.evalue_threshold,
                "query_coverage": self.coverage_threshold,
                "score": self.score_threshold,
                "identity": self.identity_threshold,
                "alignment_length": self.alignment_length_threshold}

    def analyze(self, proteins: list[Protein]) -> EvidenceAdapterResult:
        provenance = self.provenance()
        if not self.available():
            return EvidenceAdapterResult(
                self.name, "UNAVAILABLE", provenance=provenance,
                message="PHROGs/MMseqs2 unavailable or database prefix is not prepared; no PHROGs evidence was fabricated.")
        with tempfile.TemporaryDirectory(prefix="phagemine-phrogs-") as temp:
            root = Path(temp)
            fasta, query_db, result_db = root / "proteins.faa", root / "query", root / "result"
            output, tmp = root / "phrogs.tsv", root / "tmp"
            fasta.write_text("".join(f">{p.protein_id}\n{p.sequence}\n" for p in proteins))
            commands = [
                [str(self.mmseqs), "createdb", str(fasta), str(query_db)],
                [str(self.mmseqs), "search", str(query_db), str(self.database_path), str(result_db), str(tmp)],
                [str(self.mmseqs), "convertalis", str(query_db), str(self.database_path),
                 str(result_db), str(output), "--format-output", MMSEQS_FORMAT],
            ]
            provenance["search_commands"] = commands
            provenance["search_parameters"] = {"threshold_mode": "MANUAL", "output_format": MMSEQS_FORMAT}
            for command in commands:
                result = subprocess.run(command, capture_output=True, text=True, check=False)
                if result.returncode:
                    return EvidenceAdapterResult(
                        self.name, "UNAVAILABLE", provenance=provenance,
                        message=f"MMseqs2 {command[1]} failed with exit code {result.returncode}: {(result.stderr or result.stdout).strip()}")
            evidence = self.parse_tabular(output.read_text(), proteins, provenance)
        self._mark_conflicts(evidence)
        return EvidenceAdapterResult(self.name, "REAL", evidence=evidence, provenance=provenance)

    def _classify(self, evalue: float, coverage: float | None, score: float,
                  identity: float | None, alignment_length: int) -> str:
        rejected = (
            (self.evalue_threshold is not None and evalue > self.evalue_threshold) or
            (self.coverage_threshold is not None and (coverage is None or coverage < self.coverage_threshold)) or
            (self.score_threshold is not None and score < self.score_threshold) or
            (self.identity_threshold is not None and (identity is None or identity < self.identity_threshold)) or
            (self.alignment_length_threshold is not None and alignment_length < self.alignment_length_threshold)
        )
        if rejected:
            return "REJECTED"
        configured = any(value is not None for value in (
            self.evalue_threshold, self.coverage_threshold, self.score_threshold,
            self.identity_threshold, self.alignment_length_threshold))
        return "STRONG" if configured else "MODERATE"

    @staticmethod
    def _annotation_keys(phrog_id: str) -> tuple[str, ...]:
        match = re.search(r"(\d+)$", phrog_id)
        return (phrog_id, match.group(1), f"phrog_{match.group(1)}", f"PHROG{match.group(1)}") if match else (phrog_id,)

    def _annotations(self) -> dict[str, dict[str, Any]]:
        if not self.annotations_path or not self.annotations_path.exists():
            return {}
        opener = gzip.open if self.annotations_path.suffix == ".gz" else open
        with opener(self.annotations_path, "rt", encoding="utf-8", errors="replace") as handle:
            result: dict[str, dict[str, Any]] = {}
            for row in csv.DictReader(handle, delimiter="\t"):
                phrog_id = row.get("phrog") or row.get("phrog_id") or row.get("PHROG") or row.get("PHROG_ID") or row.get("GroupName") or row.get("id")
                if phrog_id:
                    for key in self._annotation_keys(phrog_id):
                        result[key] = row
            return result

    def parse_tabular(self, text: str, proteins: list[Protein],
                      provenance: dict[str, Any] | None = None) -> list[Evidence]:
        protein_ids = {p.protein_id for p in proteins}
        annotations = self._annotations()
        records: list[Evidence] = []
        for line in text.splitlines():
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != len(MMSEQS_FIELDS):
                continue
            row = dict(zip(MMSEQS_FIELDS, fields))
            try:
                protein_id, phrog_id = row["query"], row["target"]
                identity, alignment_length = float(row["fident"]), int(row["alnlen"])
                qstart, qend, qlen = int(row["qstart"]), int(row["qend"]), int(row["qlen"])
                tstart, tend, tlen = int(row["tstart"]), int(row["tend"]), int(row["tlen"])
                evalue, score = float(row["evalue"]), float(row["bits"])
                qcov = float(row["qcov"]) if row["qcov"] else None
                tcov = float(row["tcov"]) if row["tcov"] else None
            except (ValueError, KeyError):
                continue
            if protein_id not in protein_ids:
                continue
            annotation = next((annotations[key] for key in self._annotation_keys(phrog_id) if key in annotations), {})
            category = annotation.get("category") or annotation.get("functional_category") or annotation.get("FunctionalCategory")
            description = annotation.get("annot") or annotation.get("annotation") or annotation.get("description") or annotation.get("function")
            unknown = not description or description.strip().lower() in {"unknown", "unknown function", "hypothetical protein", "na", "n/a"}
            strength = self._classify(evalue, qcov, score, identity, alignment_length)
            metrics = {
                "phrog_id": phrog_id, "query_protein_id": protein_id,
                "annotation": description, "functional_category": category,
                # Stable names plus legacy aliases retained for persisted/client compatibility.
                "percent_identity": identity, "bit_score": score,
                "mmseqs_score": score, "evalue": evalue, "sequence_identity": identity,
                "alignment_length": alignment_length, "query_coverage": qcov,
                "target_coverage": tcov, "profile_coverage": tcov,
                "query_coordinates": {"start": qstart, "end": qend},
                "target_coordinates": {"start": tstart, "end": tend},
                "query_length": qlen, "target_length": tlen,
                "raw_mmseqs_row": line, "unknown_function": unknown, "conflict": False,
            }
            supports = strength != "REJECTED"
            records.append(Evidence(
                "phage_orthology",
                f"PHROGs phage orthology evidence for {phrog_id}; evidence only, not an automatic product assignment.",
                EvidenceLevel.COMPUTATIONAL if supports else EvidenceLevel.WEAK,
                "PHROGs", self.database_version, status="REAL", supports=supports,
                metrics=metrics, identifier=phrog_id,
                threshold={**self._thresholds(), "threshold_mode": "MANUAL"},
                coordinates={"start": qstart, "end": qend},
                provenance={**(provenance or {}), "protein_id": protein_id,
                            "phrogs_version": self.database_version, "threshold_mode": "MANUAL"},
                evidence_strength=strength, family_name=phrog_id,
                # Unknown PHROGs support orthology but never fabricate a function.
                description=None if unknown else description))
        return records

    @staticmethod
    def _mark_conflicts(records: list[Evidence]) -> None:
        by_protein: dict[str, set[str]] = {}
        for record in records:
            if record.supports and record.evidence_strength == "STRONG" and record.description:
                by_protein.setdefault(record.provenance.get("protein_id", ""), set()).add(record.description)
        for record in records:
            if len(by_protein.get(record.provenance.get("protein_id", ""), set())) > 1:
                record.metrics["conflict"] = True
                record.provenance["conflict"] = True


# Import compatibility only; execution and evidence provenance are MMseqs2.
PHROGSHMMAdapter = PHROGSMMseqsAdapter
