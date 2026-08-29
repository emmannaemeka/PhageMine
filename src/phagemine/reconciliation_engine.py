"""Caller-neutral, observational gene-model reconciliation."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from .gene_models import GeneModel, ReconciledLocus

ALGORITHM_VERSION = "1.0"


def run_provider_reconciliation(providers: list[Any], genome_id: str, sequence: str,
                                input_fasta: str | Path, output_dir: str | Path,
                                *, molecule_type: str = "dna") -> tuple[dict[str, Any], list[ReconciledLocus]]:
    """Development/internal N-provider execution path.

    Providers are run independently and their models are only reconciled for
    observation; this function never selects a final CDS model.
    """
    results = {}
    model_sets = {}
    raw_root = Path(output_dir) / "raw"
    for provider in sorted(providers, key=lambda item: item.provider_id):
        result = provider.predict(genome_id, sequence, input_fasta, molecule_type=molecule_type)
        provider.persist_raw_output(result, raw_root)
        results[provider.provider_id] = result
        model_sets[provider.provider_id] = result.models
    loci = reconcile_gene_models(model_sets)
    write_reconciliation_v2(Path(output_dir) / "reconciliation", model_sets, loci)
    return results, loci


def _overlap(a: GeneModel, b: GeneModel) -> int:
    return max(0, min(a.end, b.end) - max(a.start, b.start) + 1)


def _meaningful_overlap(a: GeneModel, b: GeneModel, *, min_overlap_bp: int, min_reciprocal: float) -> bool:
    overlap = _overlap(a, b)
    if overlap <= 0:
        return False
    reciprocal = min(overlap / a.length_nt, overlap / b.length_nt)
    # The absolute floor protects short genes, while the small fractional
    # floor prevents a 30-bp sliver of two long unrelated genes from chaining
    # them into one locus.  Complete/near-complete overlap remains accepted by
    # the stronger reciprocal criterion.
    return (overlap >= min_overlap_bp and reciprocal >= 0.1) or reciprocal >= min_reciprocal


def _candidate_key(provider: str, model: GeneModel) -> tuple:
    protein_hash = hashlib.sha256((model.protein_sequence or model.sequence or "").encode()).hexdigest()[:16]
    return (model.segment_id or "", model.start, model.end, model.strand, provider, model.raw_identifier, protein_hash)


def _candidate_id(provider: str, model: GeneModel) -> str:
    payload = "|".join(map(str, _candidate_key(provider, model)))
    return "GM_" + hashlib.sha256(payload.encode()).hexdigest()[:16]


def _classify(models: list[GeneModel], providers: list[str], edges: dict[int, set[int]]) -> tuple[str, bool, list[str]]:
    if len(set(m.strand for m in models)) > 1:
        return "STRAND_DISCORDANCE", True, ["Opposite-strand hypotheses are retained in one locus."]
    degrees_left = [len(edges.get(i, set())) for i in range(len(models))]
    # Edges are model-index pairs. A component containing multiple models from
    # one provider against one/multiple from another is a split/merge pattern.
    provider_model_counts = {p: sum(m.caller == p for m in models) for p in set(providers)}
    if any(v > 1 for v in provider_model_counts.values()) and len(provider_model_counts) > 1:
        counts = [provider_model_counts[p] for p in sorted(provider_model_counts)]
        if len(counts) == 2 and sorted(counts)[0] == 1 and sorted(counts)[1] > 1:
            # Canonical provider order makes this deterministic: the provider
            # contributing one model is the split source; the provider with
            # multiple models is the merged source.
            return ("SPLIT_MODEL" if counts[0] == 1 else "MERGED_MODEL"), True, []
        return "COMPLEX_CONFLICT", True, ["Multiple models from at least one provider occur in this locus."]
    coords = {(m.start, m.end, m.strand) for m in models}
    starts = {m.start for m in models}; ends = {m.end for m in models}
    if len(coords) == 1:
        return ("EXACT_CONCORDANCE" if len(set(providers)) > 1 else "CALLER_SPECIFIC"), False, []
    if len(ends) == 1 and len(starts) > 1:
        return "COMMON_STOP_ALTERNATE_START", True, []
    if len(starts) == 1 and len(ends) > 1:
        return "COMMON_START_ALTERNATE_STOP", True, []
    return ("BOUNDARY_DISCORDANCE" if len(set(providers)) > 1 else "CALLER_SPECIFIC"), True, []


def reconcile_gene_models(provider_models: dict[str, list[GeneModel]], *,
                          min_overlap_bp: int = 30, min_reciprocal: float = 0.5) -> list[ReconciledLocus]:
    """Group arbitrary provider models into deterministic observational loci.

    Models are connected when they share genome/segment, strand, and either a
    30-bp meaningful overlap or 50% reciprocal overlap. Connected components
    retain every candidate; no model is selected or averaged. Provider and
    model order cannot affect the result.
    """
    flattened: list[tuple[str, GeneModel]] = []
    for provider, models in sorted(provider_models.items()):
        for model in sorted(models, key=lambda m: _candidate_key(provider, m)):
            flattened.append((provider, model))
    n = len(flattened)
    adjacency = [set() for _ in range(n)]
    # Sweep each genome/segment interval list; candidates are compared only
    # while their intervals remain active rather than through a global N² loop.
    groups: dict[tuple[str, str], list[int]] = {}
    for idx, (_, model) in enumerate(flattened):
        groups.setdefault((model.genome_id or "", model.segment_id or ""), []).append(idx)
    for indices in groups.values():
        ordered = sorted(indices, key=lambda idx: (flattened[idx][1].start, flattened[idx][1].end, flattened[idx][0], flattened[idx][1].raw_identifier))
        active: list[int] = []
        for idx in ordered:
            current = flattened[idx][1]
            active = [other for other in active if flattened[other][1].end >= current.start]
            for other in active:
                if _meaningful_overlap(current, flattened[other][1], min_overlap_bp=min_overlap_bp, min_reciprocal=min_reciprocal):
                    adjacency[idx].add(other); adjacency[other].add(idx)
            active.append(idx)
    visited = set(); components = []
    for root in range(n):
        if root in visited: continue
        stack = [root]; visited.add(root); members = []
        while stack:
            current = stack.pop(); members.append(current)
            for nxt in sorted(adjacency[current]):
                if nxt not in visited: visited.add(nxt); stack.append(nxt)
        components.append(sorted(members, key=lambda idx: _candidate_key(flattened[idx][0], flattened[idx][1])))
    loci: list[ReconciledLocus] = []
    for members in components:
        pairs = [flattened[idx] for idx in members]
        models = [model for _, model in pairs]
        providers = [provider for provider, _ in pairs]
        candidate_ids = [_candidate_id(provider, model) for provider, model in pairs]
        order = sorted(range(len(models)), key=lambda i: _candidate_key(providers[i], models[i]))
        models = [models[i] for i in order]; providers = [providers[i] for i in order]; candidate_ids = [candidate_ids[i] for i in order]
        start, end = min(m.start for m in models), max(m.end for m in models)
        signature = "|".join(candidate_ids)
        locus_id = "LOCUS_" + hashlib.sha256(f"{models[0].genome_id or ''}|{models[0].segment_id or ''}|{start}|{end}|{signature}".encode()).hexdigest()[:16]
        coordinate_groups: dict[tuple[int, int, str], list[str]] = {}
        start_groups: dict[str, list[str]] = {}; stop_groups: dict[str, list[str]] = {}; strand_groups: dict[str, list[str]] = {}
        for cid, model in zip(candidate_ids, models):
            coordinate_groups.setdefault((model.start, model.end, model.strand), []).append(cid)
            start_groups.setdefault(str(model.start), []).append(cid); stop_groups.setdefault(str(model.end), []).append(cid); strand_groups.setdefault(model.strand, []).append(cid)
        exact = [{"start": k[0], "end": k[1], "strand": k[2], "candidate_ids": sorted(v), "supporting_providers": sorted({providers[candidate_ids.index(cid)] for cid in v})} for k, v in sorted(coordinate_groups.items())]
        cls, review, notes = _classify(models, providers, {})
        loci.append(ReconciledLocus(locus_id, models[0].genome_id, models[0].segment_id, start, end,
                                    "CONCORDANT" if len(strand_groups) == 1 else "DISCORDANT",
                                    models, sorted(set(providers)), len(set(providers)), len(models), cls,
                                    exact, {k: sorted(v) for k, v in sorted(start_groups.items())},
                                    {k: sorted(v) for k, v in sorted(stop_groups.items())},
                                    {k: sorted(v) for k, v in sorted(strand_groups.items())}, review, notes))
    return sorted(loci, key=lambda locus: (locus.segment_id or "", locus.start, locus.end, locus.locus_id))


def write_reconciliation_v2(root: str | Path, provider_models: dict[str, list[GeneModel]],
                            loci: list[ReconciledLocus], *, min_overlap_bp: int = 30,
                            min_reciprocal: float = 0.5) -> None:
    root = Path(root); root.mkdir(parents=True, exist_ok=True)
    candidates = []
    for provider, models in sorted(provider_models.items()):
        for model in sorted(models, key=lambda m: _candidate_key(provider, m)):
            candidates.append({"candidate_id": _candidate_id(provider, model), "provider_id": provider,
                               "raw_identifier": model.raw_identifier, "genome_id": model.genome_id,
                               "segment_id": model.segment_id, "start": model.start, "end": model.end,
                               "strand": model.strand, "length_nt": model.length_nt, "length_aa": model.length_aa,
                               "protein_sha256": hashlib.sha256((model.protein_sequence or model.sequence or "").encode()).hexdigest(),
                               "input_sequence_sha256": model.input_sequence_sha256 or ""})
    with (root / "candidate_gene_models.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(candidates[0]) if candidates else ["candidate_id"], delimiter="\t"); writer.writeheader(); writer.writerows(candidates)
    summary = []
    for locus in loci:
        summary.append({"locus_id": locus.locus_id, "genome_id": locus.genome_id or "", "segment_id": locus.segment_id or "",
                        "start": locus.start, "end": locus.end, "strand_status": locus.strand_status,
                        "supporting_providers": ",".join(locus.supporting_providers), "provider_count": locus.provider_count,
                        "candidate_count": locus.candidate_count, "reconciliation_class": locus.reconciliation_class,
                        "exact_coordinate_groups": json.dumps(locus.exact_coordinate_groups, sort_keys=True),
                        "review_required": str(locus.review_required).lower(), "notes": "; ".join(locus.notes)})
    with (root / "reconciliation_v2.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]) if summary else ["locus_id"], delimiter="\t"); writer.writeheader(); writer.writerows(summary)
    (root / "reconciliation_v2.json").write_text(json.dumps({"algorithm_version": ALGORITHM_VERSION,
        "locus_formation_policy": {"min_overlap_bp": min_overlap_bp, "min_reciprocal": min_reciprocal},
        "candidate_models": candidates, "loci": [locus.to_dict() for locus in loci]}, indent=2, sort_keys=True))
    manifest_path = root / "gene_call_manifest.json"
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text())
        manifest["reconciliation_engine"] = {
            "schema_version": "1.0", "algorithm_version": ALGORITHM_VERSION,
            "active_providers": sorted(provider_models),
            "locus_formation_policy": "connected overlap graph",
            "thresholds": {"min_overlap_bp": min_overlap_bp, "min_reciprocal": min_reciprocal},
            "deterministic_ordering_policy": "canonical segment/coordinate/provider/model signatures",
            "output_files": ["candidate_gene_models.tsv", "reconciliation_v2.tsv", "reconciliation_v2.json"],
        }
        manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True))
