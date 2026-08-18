"""Generic interfaces for optional, provenance-preserving evidence adapters."""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any
from enum import Enum

from .models import Evidence, Protein


class EvidenceState(str, Enum):
    NOT_REQUESTED = "NOT_REQUESTED"
    READY = "READY"
    RUNNING = "RUNNING"
    SUCCESS_NO_HIT = "SUCCESS_NO_HIT"
    SUCCESS_WITH_HITS = "SUCCESS_WITH_HITS"
    UNAVAILABLE = "UNAVAILABLE"
    FAILED = "FAILED"
    STALE = "STALE"


@dataclass
class EvidenceAdapterResult:
    adapter: str
    status: str
    evidence: list[Evidence] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)
    message: str | None = None

    @property
    def state(self) -> EvidenceState:
        if self.status in {"NOT_REQUESTED", "READY", "RUNNING", "SUCCESS_NO_HIT", "SUCCESS_WITH_HITS", "UNAVAILABLE", "FAILED", "STALE"}:
            return EvidenceState(self.status)
        if self.status == "REAL":
            return EvidenceState.SUCCESS_WITH_HITS if self.evidence else EvidenceState.SUCCESS_NO_HIT
        return EvidenceState.UNAVAILABLE if self.status == "UNAVAILABLE" else EvidenceState.FAILED


class EvidenceAdapter(ABC):
    name: str

    @abstractmethod
    def version(self) -> str: ...

    @abstractmethod
    def available(self) -> bool: ...

    @abstractmethod
    def analyze(self, proteins: list[Protein]) -> EvidenceAdapterResult: ...

    @abstractmethod
    def provenance(self) -> dict[str, Any]: ...
