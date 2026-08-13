from __future__ import annotations

from .models import Evidence, EvidenceLevel, Protein
from .evidence import EvidenceAdapter, EvidenceAdapterResult


class MockEvidenceBackend(EvidenceAdapter):
    """Deterministic demonstration evidence; never use as biological database output."""
    name = "mock-phage-evidence"

    def version(self) -> str:
        return "demo-1"

    def available(self) -> bool:
        return True

    def provenance(self) -> dict:
        return {"adapter": self.name, "adapter_version": self.version(), "status": "MOCK", "message": "Demonstration evidence only; not biological evidence."}

    def analyze(self, proteins: list[Protein]) -> EvidenceAdapterResult:
        self.annotate(proteins)
        return EvidenceAdapterResult(self.name, "MOCK", provenance=self.provenance())

    def annotate(self, proteins: list[Protein]) -> None:
        for index, protein in enumerate(proteins, 1):
            # Retain substantial unknown space; evidence is deterministic for reproducible fixtures.
            selector = index % 5
            if selector == 1:
                protein.annotation = "Putative structural protein (MOCK)"
                protein.annotation_level = EvidenceLevel.CURATED
                protein.evidence.append(Evidence("similarity", "Mock curated-like homolog supports a structural-protein family.", EvidenceLevel.CURATED, self.name, self.version(), metrics={"identity": 0.64, "coverage": 0.91, "e_value": "1e-30"}))
            elif selector == 2:
                protein.annotation = "Uncharacterized protein with weak domain signal (MOCK)"
                protein.annotation_level = EvidenceLevel.WEAK
                protein.evidence.append(Evidence("domain", "Mock weak profile signal; insufficient for exact function.", EvidenceLevel.WEAK, self.name, self.version(), metrics={"coverage": 0.42, "e_value": "2e-4"}))
            elif selector == 3:
                protein.annotation = "Hypothetical protein"
                protein.annotation_level = EvidenceLevel.HYPOTHESIS
            elif selector == 4:
                protein.annotation = "Putative DNA metabolism protein (MOCK)"
                protein.annotation_level = EvidenceLevel.COMPUTATIONAL
                protein.evidence.append(Evidence("domain", "Mock conserved catalytic-family profile supports broad DNA-metabolism association.", EvidenceLevel.COMPUTATIONAL, self.name, self.version(), metrics={"coverage": 0.72, "e_value": "3e-12"}))
            else:
                protein.annotation = "Hypothetical protein"
                protein.annotation_level = EvidenceLevel.HYPOTHESIS
                protein.evidence.append(Evidence("quality", "Mock low-complexity warning reduces interpretability.", EvidenceLevel.WEAK, self.name, self.version(), supports=False, metrics={"low_complexity_fraction": 0.31}))
