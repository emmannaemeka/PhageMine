"""Generic gene-model provider contract and registry.

Providers are prediction engines only.  They do not select final models or
perform reconciliation; those policies remain in later pipeline layers.
"""
from __future__ import annotations

import hashlib
import csv
import json
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any

from .gene_models import GeneModel
from .gene_prediction import PHANOTATEPredictor
from .reconciliation import ProdigalPredictor


class MoleculeType(str, Enum):
    DNA = "dna"
    RNA = "rna"


class ProviderError(RuntimeError):
    """Base class for provider-level, machine-classifiable failures."""


class ProviderUnavailable(ProviderError):
    code = "PROVIDER_UNAVAILABLE"


class UnsupportedMoleculeType(ProviderError):
    code = "UNSUPPORTED_MOLECULE_TYPE"


class ProviderExecutionFailure(ProviderError):
    code = "EXECUTION_FAILURE"


class NoParseableCalls(ProviderError):
    code = "NO_PARSEABLE_CALLS"


class InvalidCallerOutput(ProviderError):
    code = "INVALID_CALLER_OUTPUT"


@dataclass(frozen=True)
class ProviderCapabilities:
    supported_molecule_types: tuple[str, ...]
    supports_segmented_input: bool = False
    supports_alternative_genetic_codes: bool = False
    supports_overlapping_orfs: bool = True
    supports_small_orfs: bool = True
    external_executable_required: bool = True
    native_python_provider: bool = False


@dataclass
class GenePredictionResult:
    provider_id: str
    provider_name: str
    provider_version: str
    models: list[GeneModel] = field(default_factory=list)
    raw_output_paths: list[str] = field(default_factory=list)
    command: list[str] | None = None
    parameters: dict[str, Any] = field(default_factory=dict)
    input_sequence_sha256: str = ""
    molecule_type: str = MoleculeType.DNA.value
    segment_id: str | None = None
    status: str = "SUCCESS"
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {
            "provider_id": self.provider_id,
            "provider_name": self.provider_name,
            "provider_version": self.provider_version,
            "models": [model.to_dict() for model in self.models],
            "raw_output_paths": list(self.raw_output_paths),
            "command": self.command,
            "parameters": dict(self.parameters),
            "input_sequence_sha256": self.input_sequence_sha256,
            "molecule_type": self.molecule_type,
            "segment_id": self.segment_id,
            "status": self.status,
            "warnings": list(self.warnings),
        }


class GeneModelProvider(ABC):
    """Prediction-engine interface independent of any particular algorithm."""

    provider_id: str
    name: str
    capabilities: ProviderCapabilities

    @abstractmethod
    def version(self) -> str: ...

    def available(self) -> bool:
        return True

    def supports_molecule_type(self, molecule_type: str | MoleculeType) -> bool:
        value = molecule_type.value if isinstance(molecule_type, MoleculeType) else str(molecule_type).lower()
        return value in self.capabilities.supported_molecule_types

    @abstractmethod
    def parameters(self) -> dict[str, Any]: ...

    @abstractmethod
    def predict(
        self, genome_id: str, sequence: str, input_fasta: str | Path | None = None,
        *, molecule_type: str | MoleculeType = MoleculeType.DNA,
        segment_id: str | None = None,
    ) -> GenePredictionResult: ...

    def raw_output_metadata(self) -> dict[str, Any]:
        return {}

    def persist_raw_output(self, result: GenePredictionResult, output_dir: str | Path) -> list[str]:
        """Persist provider-native output when the provider has one.

        Providers must not invent a representation that their API did not
        produce.  The default therefore records no files.
        """
        return []


def _model_from_protein(protein, *, provider_id: str, version: str, command, input_sha: str, source_file):
    params = dict(getattr(protein, "gene_call_parameters", {}) or {})
    return GeneModel(
        caller=provider_id, identifier=protein.protein_id, start=protein.start,
        end=protein.end, strand=protein.strand, sequence=protein.sequence,
        caller_version=version, command=command, options=params,
        genome_id=protein.genome_id, start_codon=protein.start_codon,
        stop_codon=protein.stop_codon, cds_sequence=protein.cds,
        protein_sequence=protein.sequence, input_sequence_sha256=input_sha,
        raw_start=params.get("raw_start", protein.start), raw_end=params.get("raw_end", protein.end),
        raw_strand=params.get("reported_strand", protein.strand), source_file=str(source_file) if source_file else None,
    )


class PHANOTATEProvider(GeneModelProvider):
    provider_id = "phanotate"
    name = "PHANOTATE"
    capabilities = ProviderCapabilities((MoleculeType.DNA.value,), supports_segmented_input=False,
                                        supports_alternative_genetic_codes=False,
                                        supports_overlapping_orfs=True, supports_small_orfs=True,
                                        external_executable_required=True, native_python_provider=False)

    def __init__(self, executable: str | None = None, extra_args: list[str] | None = None):
        self.predictor = PHANOTATEPredictor(executable, extra_args)

    def version(self) -> str:
        return self.predictor.version()

    def available(self) -> bool:
        return self.predictor.available()

    def parameters(self) -> dict[str, Any]:
        return self.predictor.parameters()

    def predict(self, genome_id, sequence, input_fasta=None, *, molecule_type=MoleculeType.DNA, segment_id=None):
        if not self.supports_molecule_type(molecule_type):
            raise UnsupportedMoleculeType(f"{self.provider_id} does not support {molecule_type}")
        if not self.available():
            raise ProviderUnavailable("PHANOTATE executable is unavailable")
        try:
            proteins = self.predictor.predict(genome_id, sequence, input_fasta)
        except ValueError as exc:
            raise InvalidCallerOutput(str(exc)) from exc
        except RuntimeError as exc:
            raise ProviderExecutionFailure(str(exc)) from exc
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        models = [_model_from_protein(p, provider_id=self.provider_id, version=self.version(), command=self.predictor.last_command, input_sha=digest, source_file=input_fasta) for p in proteins]
        if not models:
            raise NoParseableCalls("PHANOTATE produced no parseable models")
        return GenePredictionResult(self.provider_id, self.name, self.version(), models,
                                    command=self.predictor.last_command, parameters=self.parameters(),
                                    input_sequence_sha256=digest, molecule_type="dna", segment_id=segment_id)

    def raw_output_metadata(self):
        return {"stdout": self.predictor.last_raw_output, "command": self.predictor.last_command}

    def persist_raw_output(self, result, output_dir):
        root = Path(output_dir); root.mkdir(parents=True, exist_ok=True)
        raw = root / "phanotate.raw.txt"; raw.write_text(self.predictor.last_raw_output)
        tsv = root / "phanotate.tsv"
        with tsv.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["caller", "raw_identifier", "start", "end", "strand", "length_nt", "length_aa"])
            for model in result.models:
                writer.writerow([self.provider_id, model.raw_identifier, model.start, model.end, model.strand, model.length_nt, model.length_aa])
        result.raw_output_paths = [str(raw), str(tsv)]
        return [str(raw), str(tsv)]


class ProdigalProvider(GeneModelProvider):
    provider_id = "prodigal"
    name = "Prodigal"
    capabilities = ProviderCapabilities((MoleculeType.DNA.value,), supports_segmented_input=False,
                                        supports_alternative_genetic_codes=False,
                                        supports_overlapping_orfs=True, supports_small_orfs=False,
                                        external_executable_required=True, native_python_provider=False)

    def __init__(self, executable: str | None = None):
        self.predictor = ProdigalPredictor(executable)

    def version(self) -> str:
        return self.predictor.version()

    def available(self) -> bool:
        return bool(self.predictor.executable and Path(self.predictor.executable).exists())

    def parameters(self):
        return self.predictor.parameters()

    def predict(self, genome_id, sequence, input_fasta=None, *, molecule_type=MoleculeType.DNA, segment_id=None):
        if not self.supports_molecule_type(molecule_type):
            raise UnsupportedMoleculeType(f"{self.provider_id} does not support {molecule_type}")
        if input_fasta is None:
            raise ProviderExecutionFailure("Prodigal requires an input FASTA path")
        if not self.available():
            raise ProviderUnavailable("Prodigal executable is unavailable")
        try:
            models = self.predictor.predict(input_fasta, sequence)
        except RuntimeError as exc:
            raise ProviderExecutionFailure(str(exc)) from exc
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        normalized = []
        for model in models:
            normalized.append(GeneModel(
                caller=self.provider_id, identifier=model.identifier, start=model.start, end=model.end,
                strand=model.strand, caller_version=self.version(), command=self.predictor.last_command,
                options=self.parameters(), genome_id=genome_id, segment_id=segment_id,
                input_sequence_sha256=digest, raw_start=model.start, raw_end=model.end,
                raw_strand=model.strand, source_file=str(input_fasta),
            ))
        if not normalized:
            raise NoParseableCalls("Prodigal produced no parseable models")
        return GenePredictionResult(self.provider_id, self.name, self.version(), normalized,
                                    command=self.predictor.last_command, parameters=self.parameters(),
                                    input_sequence_sha256=digest, molecule_type="dna", segment_id=segment_id)

    def raw_output_metadata(self):
        return {"gff": self.predictor.last_raw_gff, "stdout": self.predictor.last_stdout,
                "stderr": self.predictor.last_stderr, "command": self.predictor.last_command}

    def persist_raw_output(self, result, output_dir):
        root = Path(output_dir); root.mkdir(parents=True, exist_ok=True)
        gff = root / "prodigal.gff"; gff.write_text(self.predictor.last_raw_gff)
        tsv = root / "prodigal.tsv"
        with tsv.open("w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["caller", "raw_identifier", "start", "end", "strand", "length_nt"])
            for model in result.models:
                writer.writerow([self.provider_id, model.raw_identifier, model.start, model.end, model.strand, model.length_nt])
        result.raw_output_paths = [str(gff), str(tsv)]
        return [str(gff), str(tsv)]


class PyrodigalProvider(GeneModelProvider):
    """Native Pyrodigal DNA provider; predictions remain observational in Step 3."""

    provider_id = "pyrodigal"
    name = "Pyrodigal"
    capabilities = ProviderCapabilities((MoleculeType.DNA.value,), supports_segmented_input=True,
                                        supports_alternative_genetic_codes=False,
                                        supports_overlapping_orfs=True, supports_small_orfs=False,
                                        external_executable_required=False, native_python_provider=True)

    def __init__(self, *, meta: bool = True, closed: bool = False, mask: bool = False,
                 min_gene: int = 90, min_edge_gene: int = 60, max_overlap: int = 60,
                 translation_table: int = 11):
        self.meta = meta
        self.closed = closed
        self.mask = mask
        self.min_gene = min_gene
        self.min_edge_gene = min_edge_gene
        self.max_overlap = max_overlap
        self.translation_table = translation_table
        self._finder = None
        self._raw_records: list[dict[str, Any]] = []

    def _module(self):
        try:
            import pyrodigal
            return pyrodigal
        except ImportError as exc:
            raise ProviderUnavailable("Pyrodigal is required for the pyrodigal provider") from exc

    def version(self) -> str:
        try:
            module = self._module()
        except ProviderUnavailable:
            return "unavailable"
        return str(getattr(module, "__version__", "unknown"))

    def available(self) -> bool:
        try:
            self._module()
            return True
        except ProviderUnavailable:
            return False

    def parameters(self) -> dict[str, Any]:
        return {"mode": "meta" if self.meta else "single", "meta": self.meta,
                "closed": self.closed, "mask": self.mask, "min_gene": self.min_gene,
                "min_edge_gene": self.min_edge_gene, "max_overlap": self.max_overlap,
                "translation_table": "provider-native-meta" if self.meta else self.translation_table,
                "training": "metagenomic pre-trained models" if self.meta else "per-sequence training"}

    def _new_finder(self):
        module = self._module()
        try:
            return module.GeneFinder(meta=self.meta, closed=self.closed, mask=self.mask,
                                     min_gene=self.min_gene, min_edge_gene=self.min_edge_gene,
                                     max_overlap=self.max_overlap)
        except (TypeError, ValueError) as exc:
            raise InvalidCallerOutput(f"Could not configure Pyrodigal: {exc}") from exc

    def predict(self, genome_id, sequence, input_fasta=None, *, molecule_type=MoleculeType.DNA, segment_id=None):
        if not self.supports_molecule_type(molecule_type):
            raise UnsupportedMoleculeType(f"{self.provider_id} does not support {molecule_type}")
        if not sequence or any(base.upper() not in {"A", "C", "G", "T", "N"} for base in sequence):
            raise InvalidCallerOutput("Pyrodigal input must contain only A/C/G/T/N")
        finder = self._new_finder()
        if not self.meta:
            try:
                finder.train(sequence, translation_table=self.translation_table)
            except (ValueError, RuntimeError) as exc:
                raise InvalidCallerOutput(f"Pyrodigal training failed: {exc}") from exc
        try:
            genes = finder.find_genes(sequence)
        except (TypeError, RuntimeError, MemoryError) as exc:
            raise ProviderExecutionFailure(f"Pyrodigal prediction failed: {exc}") from exc
        digest = hashlib.sha256(sequence.encode()).hexdigest()
        models: list[GeneModel] = []
        self._raw_records = []
        for index, gene in enumerate(genes, 1):
            # Pyrodigal exposes zero-based inclusive begin/end coordinates.
            raw_start, raw_end = int(gene.begin), int(gene.end)
            start, end = min(raw_start, raw_end) + 1, max(raw_start, raw_end) + 1
            strand = "+" if int(gene.strand) == 1 else "-"
            cds = str(gene.sequence()).upper()
            # In meta mode Pyrodigal chooses the translation table as part of
            # its pre-trained model; preserve that effective table rather than
            # relabelling the protein with an incompatible requested table.
            protein = str(gene.translate(include_stop=False, strict=False))
            raw_id = f"PYRODIGAL_{index:06d}"
            frame = ((start - 1) % 3) + 1 if strand == "+" else ((end - 1) % 3) + 1
            record = {"raw_identifier": raw_id, "begin": raw_start, "end": raw_end,
                      "strand": int(gene.strand), "partial_begin": bool(gene.partial_begin),
                      "partial_end": bool(gene.partial_end), "start_type": str(gene.start_type),
                      "translation_table": int(gene.translation_table), "sequence": cds,
                      "protein_sequence": protein}
            self._raw_records.append(record)
            models.append(GeneModel(
                caller=self.provider_id, identifier=raw_id, start=start, end=end,
                strand=strand, sequence=protein, frame=frame, caller_version=self.version(),
                command=["pyrodigal.GeneFinder.find_genes"], options=self.parameters(),
                genome_id=genome_id, segment_id=segment_id,
                start_codon=cds[:3] if len(cds) >= 3 else None,
                stop_codon=cds[-3:] if len(cds) >= 3 else None,
                cds_sequence=cds, protein_sequence=protein, input_sequence_sha256=digest,
                raw_start=raw_start, raw_end=raw_end, raw_strand=str(gene.strand),
                coordinate_system="pyrodigal-0-based-inclusive->phagemine-1-based-inclusive",
                source_record=genome_id, source_file=str(input_fasta) if input_fasta else None,
            ))
        if not models:
            raise NoParseableCalls("Pyrodigal produced no gene models")
        return GenePredictionResult(self.provider_id, self.name, self.version(), models,
                                    command=["pyrodigal.GeneFinder.find_genes"], parameters=self.parameters(),
                                    input_sequence_sha256=digest, molecule_type="dna", segment_id=segment_id)

    def raw_output_metadata(self):
        return {"provider": "pyrodigal", "records": list(self._raw_records), "parameters": self.parameters()}

    def persist_raw_output(self, result, output_dir):
        root = Path(output_dir)
        root.mkdir(parents=True, exist_ok=True)
        tsv = root / "pyrodigal.tsv"
        with tsv.open("w", newline="") as handle:
            columns = ["raw_identifier", "begin", "end", "strand", "partial_begin", "partial_end",
                       "start_type", "translation_table", "sequence", "protein_sequence"]
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader(); writer.writerows(self._raw_records)
        native = root / "pyrodigal.json"
        native.write_text(json.dumps({"provider_id": self.provider_id, "version": self.version(),
                                      "parameters": self.parameters(), "records": self._raw_records}, indent=2, sort_keys=True))
        result.raw_output_paths = [str(tsv), str(native)]
        return [str(tsv), str(native)]


_PROVIDERS = {"phanotate": PHANOTATEProvider, "prodigal": ProdigalProvider, "pyrodigal": PyrodigalProvider}


def register_gene_model_provider(provider_id: str, provider_factory) -> None:
    key = str(provider_id).lower()
    if not key or not callable(provider_factory):
        raise ValueError("provider_id and callable provider_factory are required")
    _PROVIDERS[key] = provider_factory


def get_gene_model_provider(provider_id: str, **kwargs) -> GeneModelProvider:
    key = str(provider_id).lower()
    try:
        return _PROVIDERS[key](**kwargs)
    except KeyError as exc:
        raise ValueError(f"Unknown gene-model provider: {provider_id}") from exc


def registered_provider_ids() -> tuple[str, ...]:
    return tuple(sorted(_PROVIDERS))


def compare_gene_model_sets(primary: list[GeneModel], secondary: list[GeneModel],
                            *, primary_id: str = "primary", secondary_id: str = "secondary") -> list[dict[str, Any]]:
    """Return a deterministic descriptive pairwise comparison.

    This helper intentionally does not select a final model.  It uses the
    existing 50% reciprocal-overlap diagnostic criterion and leaves the
    later v1.2 reconciliation policy to a separate layer.
    """
    def overlap(a, b):
        return max(0, min(a.end, b.end) - max(a.start, b.start) + 1)
    rows = []
    secondary = sorted(secondary, key=lambda m: (m.start, m.end, m.raw_identifier))
    used: set[int] = set()
    for left in sorted(primary, key=lambda m: (m.start, m.end, m.raw_identifier)):
        candidates = [(idx, right, overlap(left, right)) for idx, right in enumerate(secondary) if idx not in used and overlap(left, right) > 0]
        if not candidates:
            rows.append({"locus_id": f"PAIR_{len(rows)+1:06d}", "primary_id": left.raw_identifier,
                         "secondary_id": "", "primary_start": left.start, "primary_end": left.end,
                         "primary_strand": left.strand, "secondary_start": "", "secondary_end": "",
                         "secondary_strand": "", "overlap_bp": 0, "primary_reciprocal_overlap": 0.0,
                         "secondary_reciprocal_overlap": 0.0, "start_agreement": False,
                         "stop_agreement": False, "strand_agreement": False,
                         "comparison_class": f"{primary_id.upper()}_ONLY"})
            continue
        idx, right, ov = max(candidates, key=lambda item: (item[2] / (left.length_nt or 1), -item[1].start, item[1].raw_identifier))
        rp, rs = ov / left.length_nt, ov / right.length_nt
        if min(rp, rs) < 0.5:
            rows.append({"locus_id": f"PAIR_{len(rows)+1:06d}", "primary_id": left.raw_identifier,
                         "secondary_id": "", "primary_start": left.start, "primary_end": left.end,
                         "primary_strand": left.strand, "secondary_start": "", "secondary_end": "",
                         "secondary_strand": "", "overlap_bp": ov, "primary_reciprocal_overlap": rp,
                         "secondary_reciprocal_overlap": rs, "start_agreement": False,
                         "stop_agreement": False, "strand_agreement": False,
                         "comparison_class": f"{primary_id.upper()}_ONLY"})
            continue
        used.add(idx)
        same_strand = left.strand == right.strand
        if not same_strand:
            category = "STRAND_DISCORDANCE"
        elif left.start == right.start and left.end == right.end:
            category = "EXACT_CONCORDANCE"
        elif left.end == right.end:
            category = "ALTERNATE_START"
        elif left.start == right.start:
            category = "ALTERNATE_STOP"
        else:
            category = "START_AND_STOP_DISCORDANCE"
        rows.append({"locus_id": f"PAIR_{len(rows)+1:06d}", "primary_id": left.raw_identifier,
                     "secondary_id": right.raw_identifier, "primary_start": left.start,
                     "primary_end": left.end, "primary_strand": left.strand,
                     "secondary_start": right.start, "secondary_end": right.end,
                     "secondary_strand": right.strand, "overlap_bp": ov,
                     "primary_reciprocal_overlap": rp, "secondary_reciprocal_overlap": rs,
                     "start_agreement": left.start == right.start, "stop_agreement": left.end == right.end,
                     "strand_agreement": same_strand, "comparison_class": category})
    for idx, right in enumerate(secondary):
        if idx not in used:
            rows.append({"locus_id": f"PAIR_{len(rows)+1:06d}", "primary_id": "",
                         "secondary_id": right.raw_identifier, "primary_start": "", "primary_end": "",
                         "primary_strand": "", "secondary_start": right.start, "secondary_end": right.end,
                         "secondary_strand": right.strand, "overlap_bp": 0,
                         "primary_reciprocal_overlap": 0.0, "secondary_reciprocal_overlap": 0.0,
                         "start_agreement": False, "stop_agreement": False, "strand_agreement": False,
                         "comparison_class": f"{secondary_id.upper()}_ONLY"})
    return rows


def run_observational_comparison(genome_id: str, sequence: str, input_fasta: str | Path,
                                 output_dir: str | Path, *, phanotate_executable: str | None = None,
                                 pyrodigal_kwargs: dict[str, Any] | None = None) -> tuple[GenePredictionResult, GenePredictionResult, list[dict[str, Any]]]:
    """Development-only PHANOTATE/Pyrodigal comparison; never selects final CDSs."""
    root = Path(output_dir)
    phanotate = PHANOTATEProvider(phanotate_executable)
    pyrodigal = PyrodigalProvider(**(pyrodigal_kwargs or {}))
    first = phanotate.predict(genome_id, sequence, input_fasta)
    second = pyrodigal.predict(genome_id, sequence, input_fasta)
    phanotate.persist_raw_output(first, root / "raw")
    pyrodigal.persist_raw_output(second, root / "raw")
    rows = compare_gene_model_sets(first.models, second.models, primary_id="phanotate", secondary_id="pyrodigal")
    comparison_dir = root / "comparisons"; comparison_dir.mkdir(parents=True, exist_ok=True)
    if rows:
        with (comparison_dir / "phanotate_vs_pyrodigal.tsv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader(); writer.writerows(rows)
    return first, second, rows
