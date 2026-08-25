"""Optional local PHROGs sequence- and profile-search evidence adapters."""
from __future__ import annotations

import csv
import gzip
import importlib
import importlib.util
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


def _annotation_keys(phrog_id: str) -> tuple[str, ...]:
    match = re.search(r"(\d+)$", phrog_id)
    return (phrog_id, match.group(1), f"phrog_{match.group(1)}", f"PHROG{match.group(1)}") if match else (phrog_id,)


def _load_annotations(path: Path | None) -> dict[str, dict[str, Any]]:
    if not path or not path.exists():
        return {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        result: dict[str, dict[str, Any]] = {}
        for row in csv.DictReader(handle, delimiter="\t"):
            phrog_id = (row.get("phrog") or row.get("phrog_id") or row.get("PHROG")
                        or row.get("PHROG_ID") or row.get("GroupName") or row.get("id"))
            if phrog_id:
                for key in _annotation_keys(phrog_id):
                    result[key] = row
        return result


def _annotation_fields(annotations: dict[str, dict[str, Any]], phrog_id: str) -> tuple[str | None, str | None, bool]:
    row = next((annotations[key] for key in _annotation_keys(phrog_id) if key in annotations), {})
    category = row.get("category") or row.get("functional_category") or row.get("FunctionalCategory")
    description = row.get("annot") or row.get("annotation") or row.get("description") or row.get("function")
    unknown = not description or description.strip().lower() in {
        "unknown", "unknown function", "hypothetical protein", "na", "n/a",
    }
    return category, description, unknown


def _decode(value: Any) -> str:
    return value.decode("utf-8", errors="replace") if isinstance(value, bytes) else str(value)


class PHROGSMMseqsAdapter(EvidenceAdapter):
    name = "PHROGSMMseqsAdapter"

    def __init__(self, database_path: str | Path | None = None,
                 annotations_path: str | Path | None = None,
                 mmseqs: str | None = None, database_version: str | None = None,
                 evalue_threshold: float | None = 1e-5,
                 coverage_threshold: float | None = 0.5,
                 score_threshold: float | None = None,
                 identity_threshold: float | None = None,
                 alignment_length_threshold: int | None = None, threads: int = 1):
        self.threads = max(1, int(threads))
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
                [str(self.mmseqs), "search", str(query_db), str(self.database_path), str(result_db), str(tmp), "--threads", str(self.threads)],
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
        return _annotation_keys(phrog_id)

    def _annotations(self) -> dict[str, dict[str, Any]]:
        return _load_annotations(self.annotations_path)

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
            category, description, unknown = _annotation_fields(annotations, phrog_id)
            strength = self._classify(evalue, qcov, score, identity, alignment_length)
            metrics = {
                "phrog_id": phrog_id, "query_protein_id": protein_id,
                "annotation": description, "functional_category": category,
                # Stable names plus legacy aliases retained for persisted/client compatibility.
                "percent_identity": identity, "bit_score": score,
                "mmseqs_score": score, "evalue": evalue, "sequence_identity": identity,
                "alignment_length": alignment_length, "query_coverage": qcov,
                "target_coverage": tcov, "profile_coverage": tcov,
                "search_backend": "MMseqs2", "search_backends": ["MMseqs2"],
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


class PHROGSPyHMMERAdapter(EvidenceAdapter):
    """Search protein queries against the PHROGs profile-HMM collection.

    The adapter consumes Pharokka's existing ``all_phrogs.h3m`` directly. It
    does not download, rebuild, or silently substitute a profile database.
    """

    name = "PHROGSPyHMMERAdapter"

    def __init__(self, database_path: str | Path | None = None,
                 annotations_path: str | Path | None = None,
                 database_version: str | None = None,
                 evalue_threshold: float | None = 1e-5,
                 coverage_threshold: float | None = 0.5,
                 score_threshold: float | None = None, threads: int = 1):
        self.database_path = Path(database_path).expanduser() if database_path else None
        self.annotations_path = Path(annotations_path).expanduser() if annotations_path else None
        self.database_version = database_version or "unknown"
        self.evalue_threshold = evalue_threshold
        self.coverage_threshold = coverage_threshold
        self.score_threshold = score_threshold
        self.threads = max(1, int(threads))

    @staticmethod
    def _module() -> Any | None:
        if importlib.util.find_spec("pyhmmer") is None:
            return None
        return importlib.import_module("pyhmmer")

    def available(self) -> bool:
        return bool(self.database_path and self.database_path.is_file() and self._module() is not None)

    def version(self) -> str:
        module = self._module()
        return str(getattr(module, "__version__", "unknown")) if module else "unavailable"

    def _thresholds(self) -> dict[str, Any]:
        return {"evalue": self.evalue_threshold, "query_coverage": self.coverage_threshold,
                "score": self.score_threshold}

    def provenance(self) -> dict[str, Any]:
        return {
            "adapter": self.name, "adapter_version": "0.1.0",
            "pyhmmer_version": self.version(),
            "phrogs_hmm_path": str(self.database_path) if self.database_path else None,
            "annotations_path": str(self.annotations_path) if self.annotations_path else None,
            "phrogs_version": self.database_version, "threshold_mode": "MANUAL",
            "thresholds": self._thresholds(),
            "status": "REAL" if self.available() else "UNAVAILABLE",
        }

    def _classify(self, evalue: float, coverage: float | None, score: float) -> str:
        rejected = (
            (self.evalue_threshold is not None and evalue > self.evalue_threshold) or
            (self.coverage_threshold is not None and (coverage is None or coverage < self.coverage_threshold)) or
            (self.score_threshold is not None and score < self.score_threshold)
        )
        return "REJECTED" if rejected else "STRONG"

    def analyze(self, proteins: list[Protein]) -> EvidenceAdapterResult:
        provenance = self.provenance()
        module = self._module()
        if not self.database_path or not self.database_path.is_file():
            return EvidenceAdapterResult(
                self.name, "UNAVAILABLE", provenance=provenance,
                message="PHROGs/PyHMMER profile database is unavailable; no profile evidence was fabricated.")
        if module is None:
            return EvidenceAdapterResult(
                self.name, "UNAVAILABLE", provenance=provenance,
                message="PyHMMER is not installed; no PHROGs profile evidence was fabricated.")
        if not proteins:
            return EvidenceAdapterResult(self.name, "REAL", evidence=[], provenance=provenance)
        annotations = _load_annotations(self.annotations_path)
        lengths = {protein.protein_id: len(protein.sequence) for protein in proteins}
        records: list[Evidence] = []
        try:
            with tempfile.TemporaryDirectory(prefix="phagemine-phrogs-hmm-") as temp:
                fasta = Path(temp) / "proteins.faa"
                fasta.write_text("".join(f">{p.protein_id}\n{p.sequence}\n" for p in proteins))
                alphabet = module.easel.Alphabet.amino()
                with module.plan7.HMMFile(str(self.database_path), alphabet=alphabet) as hmms:
                    with module.easel.SequenceFile(str(fasta), digital=True, alphabet=alphabet) as sequences:
                        options = {"cpus": self.threads}
                        if self.evalue_threshold is not None:
                            options["E"] = float(self.evalue_threshold)
                        for hits in module.hmmer.hmmscan(sequences, hmms, **options):
                            protein_id = _decode(hits.query.name)
                            if protein_id not in lengths:
                                continue
                            for hit in hits:
                                if not getattr(hit, "reported", True):
                                    continue
                                phrog_id = _decode(hit.name)
                                score, evalue = float(hit.score), float(hit.evalue)
                                domain = getattr(hit, "best_domain", None)
                                alignment = getattr(domain, "alignment", None) if domain is not None else None
                                query_start = int(getattr(alignment, "target_from", 0) or 0)
                                query_end = int(getattr(alignment, "target_to", 0) or 0)
                                profile_start = int(getattr(alignment, "hmm_from", 0) or 0)
                                profile_end = int(getattr(alignment, "hmm_to", 0) or 0)
                                profile_length = int(getattr(alignment, "hmm_length", 0) or 0)
                                aligned = max(0, query_end - query_start + 1) if query_start and query_end else 0
                                coverage = aligned / lengths[protein_id] if aligned else None
                                profile_coverage = ((max(0, profile_end - profile_start + 1) / profile_length)
                                                    if profile_start and profile_end and profile_length else None)
                                category, description, unknown = _annotation_fields(annotations, phrog_id)
                                strength = self._classify(evalue, coverage, score)
                                supports = strength != "REJECTED"
                                metrics = {
                                    "phrog_id": phrog_id, "query_protein_id": protein_id,
                                    "annotation": description, "functional_category": category,
                                    "bit_score": score, "evalue": evalue,
                                    "domain_i_evalue": (float(domain.i_evalue) if domain is not None else None),
                                    "alignment_length": aligned or None, "query_coverage": coverage,
                                    "profile_coverage": profile_coverage,
                                    "query_coordinates": {"start": query_start, "end": query_end},
                                    "profile_coordinates": {"start": profile_start, "end": profile_end},
                                    "query_length": lengths[protein_id], "profile_length": profile_length or None,
                                    "unknown_function": unknown, "conflict": False,
                                    "search_backend": "PyHMMER", "search_backends": ["PyHMMER"],
                                }
                                records.append(Evidence(
                                    "phage_profile_hmm",
                                    f"PHROGs profile-HMM evidence for {phrog_id}; evidence only, not an automatic product assignment.",
                                    EvidenceLevel.COMPUTATIONAL if supports else EvidenceLevel.WEAK,
                                    "PHROGs", self.database_version, status="REAL", supports=supports,
                                    metrics=metrics, identifier=phrog_id,
                                    threshold={**self._thresholds(), "threshold_mode": "MANUAL"},
                                    coordinates={"start": query_start, "end": query_end} if query_start else None,
                                    provenance={**provenance, "protein_id": protein_id,
                                                "phrogs_version": self.database_version,
                                                "search_backend": "PyHMMER", "threshold_mode": "MANUAL"},
                                    evidence_strength=strength, family_name=phrog_id,
                                    description=None if unknown else description))
        except Exception as exc:
            provenance["error"] = f"{type(exc).__name__}: {exc}"
            return EvidenceAdapterResult(
                self.name, "UNAVAILABLE", provenance=provenance,
                message=f"PHROGs/PyHMMER failed: {type(exc).__name__}: {exc}")
        PHROGSMMseqsAdapter._mark_conflicts(records)
        return EvidenceAdapterResult(self.name, "REAL", evidence=records, provenance=provenance)


def merge_phrogs_evidence(*groups: list[Evidence]) -> list[Evidence]:
    """Merge duplicate backend hits while preserving backend corroboration."""
    merged: dict[tuple[str, str], Evidence] = {}
    for record in (item for group in groups for item in group):
        protein_id = str(record.provenance.get("protein_id") or record.metrics.get("query_protein_id") or "")
        key = (protein_id, str(record.identifier or record.family_name or ""))
        previous = merged.get(key)
        if previous is None:
            merged[key] = record
            continue
        backends = sorted(set(previous.metrics.get("search_backends", [])) |
                          set(record.metrics.get("search_backends", [])))
        previous.metrics["search_backends"] = backends
        previous.provenance["search_backends"] = backends
        previous.metrics.setdefault("corroborating_backend_hits", []).append({
            "backend": record.metrics.get("search_backend"),
            "bit_score": record.metrics.get("bit_score"),
            "evalue": record.metrics.get("evalue"),
            "query_coverage": record.metrics.get("query_coverage"),
        })
        if not previous.supports and record.supports:
            merged[key] = record
            record.metrics["search_backends"] = backends
            record.provenance["search_backends"] = backends
    result = list(merged.values())
    PHROGSMMseqsAdapter._mark_conflicts(result)
    return result


# Historical import compatibility. New profile-HMM execution uses the explicit
# PHROGSPyHMMERAdapter rather than silently aliasing MMseqs2.
PHROGSHMMAdapter = PHROGSPyHMMERAdapter
