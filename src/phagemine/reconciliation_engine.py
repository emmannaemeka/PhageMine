"""Caller-neutral, observational gene-model reconciliation."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

from .gene_models import GeneModel, ReconciledLocus
from .gene_callers import GenePredictionResult, ProviderError

ALGORITHM_VERSION = "1.0"

PRODUCTION_RELATIONSHIP_CLASSES = {
    "EXACT_CONCORDANCE": "EXACT_MATCH",
    "COMMON_STOP_ALTERNATE_START": "SAME_STOP_DIFFERENT_START",
    "COMMON_START_ALTERNATE_STOP": "NEAR_BOUNDARY_MATCH",
    "BOUNDARY_DISCORDANCE": "NEAR_BOUNDARY_MATCH",
    "CALLER_SPECIFIC": "PHANOTATE_ONLY",
    "STRAND_DISCORDANCE": "CONFLICTING_ORF",
    "COMPLEX_CONFLICT": "CONFLICTING_ORF",
    "SPLIT_MODEL": "OVERLAPPING_ALTERNATIVE",
    "MERGED_MODEL": "OVERLAPPING_ALTERNATIVE",
}


def _production_class(locus: ReconciledLocus) -> str:
    """Map the caller-neutral topology to the production vocabulary.

    This is a reporting translation only.  The generic reconciliation graph
    remains unchanged so historical consumers retain their schema.
    """
    if locus.reconciliation_class != "CALLER_SPECIFIC":
        return PRODUCTION_RELATIONSHIP_CLASSES.get(locus.reconciliation_class, "CONFLICTING_ORF")
    return "PYRODIGAL_GV_ONLY" if "phanotate" not in {p.lower() for p in locus.supporting_providers} else "PHANOTATE_ONLY"


def _production_confidence(relation: str, *, secondary_status: str) -> tuple[str, str]:
    if secondary_status != "SUCCESS":
        return "UNRESOLVED", "Pyrodigal-gv observational evidence was unavailable; reconciliation-dependent confidence is unresolved."
    if relation == "EXACT_MATCH":
        return "MODERATE", "PHANOTATE and Pyrodigal-gv agree exactly; caller agreement alone is not sufficient for HIGH confidence."
    if relation in {"SAME_STOP_DIFFERENT_START", "NEAR_BOUNDARY_MATCH"}:
        return "LOW", "Callers overlap but disagree at a boundary; the PHANOTATE model remains authoritative pending independent evidence."
    if relation == "PHANOTATE_ONLY":
        return "LOW", "The PHANOTATE model has no matching Pyrodigal-gv observation; independent biological support must be assessed separately."
    if relation == "PYRODIGAL_GV_ONLY":
        return "UNRESOLVED", "Only the observational caller proposed this locus; it is not added to the final CDS model."
    return "LOW", "The callers present a structural conflict; neither disagreement nor overlap establishes biological truth."


def write_production_observational_outputs(root: str | Path, provider_models: dict[str, list[GeneModel]],
                                           loci: list[ReconciledLocus], provider_results: dict[str, GenePredictionResult],
                                           *, final_provider: str = "phanotate") -> list[dict[str, Any]]:
    """Write the automatic PHANOTATE/Pyrodigal-gv observation layer.

    The returned records are keyed to PHANOTATE models for downstream review
    flags.  No candidate selection occurs here and PHANOTATE coordinates are
    never rewritten.
    """
    root = Path(root); root.mkdir(parents=True, exist_ok=True)
    secondary = provider_results.get("prodigal_gv")
    secondary_status = secondary.status if secondary else "FAILED_PROVIDER"
    records = []
    rows = []
    for locus in loci:
        relation = _production_class(locus)
        confidence, rationale = _production_confidence(relation, secondary_status=secondary_status)
        phanotate = next((model for model in locus.candidate_models if model.caller.lower() == final_provider), None)
        if phanotate is not None:
            records.append({"protein_id": phanotate.raw_identifier, "locus_id": locus.locus_id, "gene_call_confidence": confidence,
                            "review_flag": "REVIEW_REQUIRED" if confidence in {"LOW", "UNRESOLVED"} else "NONE",
                            "structural_CDS_confidence": confidence, "structural_relationship": relation,
                            "review_reason": rationale})
        rows.append({"locus_id": locus.locus_id, "start": locus.start, "end": locus.end,
                     "strand": locus.strand_status, "relationship": relation,
                     "phanotate_ids": ",".join(m.raw_identifier for m in locus.candidate_models if m.caller.lower() == final_provider),
                     "pyrodigal_gv_ids": ",".join(m.raw_identifier for m in locus.candidate_models if m.caller.lower() == "prodigal_gv"),
                     "structural_CDS_confidence": confidence, "rationale": rationale,
                     "independent_biological_evidence": "NOT_ASSESSED_AT_STRUCTURAL_STAGE",
                     "independent_evidence_sources": "", "independent_evidence_count": 0, "evidence_conflict": False,
                     "review_required": "true" if confidence in {"LOW", "UNRESOLVED"} else str(locus.review_required).lower()})
    fields = ["locus_id", "start", "end", "strand", "relationship", "phanotate_ids", "pyrodigal_gv_ids",
              "structural_CDS_confidence", "rationale", "independent_biological_evidence",
              "independent_evidence_sources", "independent_evidence_count", "evidence_conflict", "review_required"]
    with (root / "structural_reconciliation.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t"); writer.writeheader(); writer.writerows(rows)
    (root / "structural_reconciliation.json").write_text(json.dumps({
        "schema_version": "production-observational-v1", "final_model_policy": "PHANOTATE_AUTHORITATIVE",
        "caller_roles": {"phanotate": "PRIMARY", "prodigal_gv": "SECONDARY_OBSERVATIONAL"},
        "secondary_status": secondary_status, "relationships": rows,
    }, indent=2, sort_keys=True))
    return records


def run_automatic_observational_reconciliation(primary_models: list[GeneModel], genome_id: str, sequence: str,
                                               input_fasta: str | Path, output_dir: str | Path,
                                               *, pyrodigal_gv_kwargs: dict[str, Any] | None = None) -> tuple[dict[str, GenePredictionResult], list[ReconciledLocus], list[dict[str, Any]]]:
    """Run the production secondary observation without changing final CDSs."""
    from .gene_callers import ProdigalGVProvider
    root = Path(output_dir); raw_root = root / "raw"; raw_root.mkdir(parents=True, exist_ok=True)
    provider = ProdigalGVProvider(**(pyrodigal_gv_kwargs or {}))
    digest = hashlib.sha256(sequence.encode()).hexdigest()
    results: dict[str, GenePredictionResult] = {}
    try:
        result = provider.predict(genome_id, sequence, input_fasta, molecule_type="dna")
        provider.persist_raw_output(result, raw_root)
    except Exception as exc:  # observational failures are deliberately non-fatal
        result = GenePredictionResult(provider_id=provider.provider_id, provider_name=provider.name,
                                      provider_version=provider.version(), input_sequence_sha256=digest,
                                      molecule_type="dna", status="FAILED_PROVIDER", warnings=[str(exc)],
                                      parameters=provider.parameters())
    results[provider.provider_id] = result
    model_sets = {"phanotate": list(primary_models), "prodigal_gv": list(result.models)}
    loci = reconcile_gene_models(model_sets)
    write_reconciliation_v2(root / "reconciliation", model_sets, loci)
    records = write_production_observational_outputs(root, model_sets, loci, results)
    return results, loci, records


def finalize_production_structural_confidence(root: str | Path, records: list[dict[str, Any]], proteins: list[Any]) -> list[dict[str, Any]]:
    """Integrate independent protein evidence without changing CDS coordinates."""
    by_id = {str(p.protein_id): p for p in proteins}
    for record in records:
        protein = by_id.get(str(record.get("protein_id")))
        evidence = list(getattr(protein, "evidence", []) or []) if protein else []
        supported = [item for item in evidence if getattr(item, "supports", False)]
        strong = [item for item in supported if getattr(item, "evidence_strength", "").upper() in {"STRONG", "EXPERIMENTAL"}]
        conflicts = [item for item in evidence if not getattr(item, "supports", True) or getattr(item, "conflict", False)]
        relation = record.get("structural_relationship")
        if strong and relation == "EXACT_MATCH":
            confidence = "HIGH"
        elif supported and not conflicts and relation in {"EXACT_MATCH", "PHANOTATE_ONLY"}:
            confidence = "MODERATE"
        elif conflicts and not supported:
            confidence = "LOW"
        else:
            confidence = record.get("structural_CDS_confidence", "UNRESOLVED")
        record["structural_CDS_confidence"] = confidence
        record["gene_call_confidence"] = confidence
        record["independent_evidence_sources"] = ";".join(sorted({str(item.source) for item in evidence if getattr(item, "source", None)}))
        record["independent_evidence_count"] = len(evidence)
        record["evidence_conflict"] = bool(conflicts)
        record["review_reason"] = (record.get("review_reason", "") +
                                    " Independent biological evidence was incorporated after evidence fusion.").strip()
    path = Path(root)
    json_path = path / "structural_reconciliation.json"
    if json_path.is_file():
        payload = json.loads(json_path.read_text())
        by_locus = {row["locus_id"]: row for row in payload.get("relationships", [])}
        for record in records:
            row = by_locus.get(record.get("locus_id"))
            if row:
                row.update({key: record[key] for key in ("structural_CDS_confidence", "independent_evidence_sources", "independent_evidence_count", "evidence_conflict") if key in record})
        json_path.write_text(json.dumps(payload, indent=2, sort_keys=True))
    tsv_path = path / "structural_reconciliation.tsv"
    if tsv_path.is_file():
        with tsv_path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        by_locus = {record.get("locus_id"): record for record in records}
        for row in rows:
            record = by_locus.get(row.get("locus_id"))
            if record:
                for key in ("structural_CDS_confidence", "independent_evidence_sources", "independent_evidence_count", "evidence_conflict"):
                    row[key] = str(record.get(key, ""))
        with tsv_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]) if rows else ["locus_id"], delimiter="\t")
            writer.writeheader(); writer.writerows(rows)
    return records


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
        try:
            result = provider.predict(genome_id, sequence, input_fasta, molecule_type=molecule_type)
            provider.persist_raw_output(result, raw_root)
            model_sets[provider.provider_id] = result.models
        except ProviderError as exc:
            # Alternative callers are diagnostic.  A provider-specific input
            # incompatibility is recorded and does not suppress the primary
            # PHANOTATE stream or masquerade as consensus support.
            result = GenePredictionResult(
                provider_id=provider.provider_id, provider_name=provider.name,
                provider_version=provider.version(), molecule_type=molecule_type,
                segment_id=None, status="SKIPPED_INCOMPATIBLE" if getattr(exc, "code", "") == "INVALID_CALLER_OUTPUT" else "FAILED_PROVIDER",
                warnings=[str(exc)])
            model_sets[provider.provider_id] = []
        results[provider.provider_id] = result
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
        # Non-segmented providers may omit segment_id while another provider
        # uses the sequence/genome identifier. Canonicalize that representation
        # so equivalent calls cannot be split into duplicate loci.
        segment_key = model.segment_id or (model.genome_id if model.genome_id else "")
        groups.setdefault((model.genome_id or "", segment_key), []).append(idx)
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
        exact = []
        for k, values in sorted(coordinate_groups.items()):
            supporting = sorted({providers[candidate_ids.index(cid)] for cid in values})
            exact.append({"start": k[0], "end": k[1], "strand": k[2], "candidate_ids": sorted(values),
                          "supporting_providers": supporting,
                          "supporting_method_families": sorted({models[candidate_ids.index(cid)].method_family or providers[candidate_ids.index(cid)] for cid in values}),
                          "supporting_method_lineages": sorted({models[candidate_ids.index(cid)].method_lineage or providers[candidate_ids.index(cid)] for cid in values})})
        cls, review, notes = _classify(models, providers, {})
        families = sorted({model.method_family or provider for provider, model in zip(providers, models)})
        lineages = sorted({model.method_lineage or provider for provider, model in zip(providers, models)})
        loci.append(ReconciledLocus(
            locus_id=locus_id, genome_id=models[0].genome_id, segment_id=models[0].segment_id,
            start=start, end=end, strand_status="CONCORDANT" if len(strand_groups) == 1 else "DISCORDANT",
            candidate_models=models, supporting_providers=sorted(set(providers)),
            provider_count=len(set(providers)), method_family_count=len(families),
            method_lineage_count=len(lineages), supporting_method_families=families,
            supporting_method_lineages=lineages, candidate_count=len(models),
            reconciliation_class=cls, exact_coordinate_groups=exact,
            start_groups={k: sorted(v) for k, v in sorted(start_groups.items())},
            stop_groups={k: sorted(v) for k, v in sorted(stop_groups.items())},
            strand_groups={k: sorted(v) for k, v in sorted(strand_groups.items())},
            review_required=review, notes=notes))
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
