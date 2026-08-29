"""Generic gene-model provider contract and registry.

Providers are prediction engines only.  They do not select final models or
perform reconciliation; those policies remain in later pipeline layers.
"""
from __future__ import annotations

import hashlib
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


_PROVIDERS = {"phanotate": PHANOTATEProvider, "prodigal": ProdigalProvider}


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
