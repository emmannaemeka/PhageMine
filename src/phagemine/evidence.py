"""Generic interfaces for optional, provenance-preserving evidence adapters."""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any

from .models import Evidence, Protein


@dataclass
class EvidenceAdapterResult:
    adapter: str
    status: str
    evidence: list[Evidence] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)
    message: str | None = None


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
