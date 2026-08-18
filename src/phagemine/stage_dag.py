"""Central dependency DAG and checksum-based stage invalidation."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class Stage:
    identifier: str
    dependencies: tuple[str, ...] = ()
    resource_bound: bool = False


STAGES: tuple[Stage, ...] = (
    Stage("GENOME_VALIDATION"),
    Stage("GENE_PREDICTION", ("GENOME_VALIDATION",)),
    Stage("PFAM", ("GENE_PREDICTION",), True),
    Stage("VOGDB", ("GENE_PREDICTION",), True),
    Stage("SWISSPROT", ("GENE_PREDICTION",), True),
    Stage("PHROGS", ("GENE_PREDICTION",), True),
    Stage("EVIDENCE_FUSION", ("PFAM", "VOGDB", "SWISSPROT", "PHROGS")),
    Stage("RANKING_MINING", ("EVIDENCE_FUSION",)),
    Stage("QC_REPORTING", ("RANKING_MINING",)),
    Stage("GENBANK", ("QC_REPORTING",)),
)
STAGE_BY_ID = {stage.identifier: stage for stage in STAGES}


def downstream(stage: str) -> set[str]:
    affected = {stage}
    changed = True
    while changed:
        changed = False
        for candidate in STAGES:
            if candidate.identifier not in affected and any(dep in affected for dep in candidate.dependencies):
                affected.add(candidate.identifier)
                changed = True
    return affected


def stale_stages(previous: dict[str, Any], current_inputs: dict[str, Any]) -> list[str]:
    """Return stages invalidated by checksum/provenance changes."""
    stale: set[str] = set()
    old_inputs = previous.get("inputs", {})
    for stage in STAGES:
        old = old_inputs.get(stage.identifier, {})
        new = current_inputs.get(stage.identifier, {})
        if old.get("input_checksums") != new.get("input_checksums") or old.get("resource_checksums") != new.get("resource_checksums") or old.get("executable_version") != new.get("executable_version"):
            stale.update(downstream(stage.identifier))
    return [stage.identifier for stage in STAGES if stage.identifier in stale]


def stage_record(identifier: str, *, status: str, input_checksums=None, output_checksums=None,
                 resource_checksums=None, executable_version=None, failure_reason=None,
                 stale: bool = False, started_at=None, ended_at=None) -> dict[str, Any]:
    stage = STAGE_BY_ID[identifier]
    return {"stage_identifier": identifier, "status": status,
            "input_checksums": input_checksums or {}, "output_checksums": output_checksums or {},
            "resource_checksums": resource_checksums or {}, "executable_version": executable_version,
            "dependencies": list(stage.dependencies), "stale": stale,
            "failure_reason": failure_reason, "started_at": started_at, "ended_at": ended_at}
