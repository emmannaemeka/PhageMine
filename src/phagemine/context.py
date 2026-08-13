"""Conservative local genomic context and functional-module inference."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

MODULE_RULES_VERSION = "7B.1"
MODULE_KEYWORDS = {
    "head_and_packaging": ("capsid", "portal", "terminase", "head", "packaging"),
    "connector": ("connector", "neck"),
    "tail": ("tail", "baseplate", "fiber", "spike"),
    "lysis": ("holin", "endolysin", "lysis", "lysin"),
    "DNA_RNA_nucleotide_metabolism": ("polymerase", "helicase", "nuclease", "replication", "transcription", "rna", "dna"),
    "transcription_regulation": ("transcriptional", "repressor", "regulator", "promoter"),
    "integration_excision": ("integrase", "recombinase", "excision"),
    "host_takeover_AMG": ("toxin", "host takeover", "auxiliary metabolic", "moron"),
}


def _signal(result: dict[str, Any]) -> str | None:
    text = " ".join(str(result.get(key) or "") for key in ("proposed_function", "functional_category")).lower()
    matches = [module for module, words in MODULE_KEYWORDS.items() if any(word in text for word in words)]
    return matches[0] if len(matches) == 1 else None


def build_context(proteins, classifications: list[dict[str, Any]], window: int = 2) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    results = {item["protein_id"]: item for item in classifications}
    ordered = sorted(proteins, key=lambda p: (p.start, p.end, p.protein_id))
    signals = [_signal(results.get(p.protein_id, {})) for p in ordered]
    segments: list[dict[str, Any]] = []
    for module in sorted({s for s in signals if s}):
        anchors = [i for i, signal in enumerate(signals) if signal == module]
        groups: list[list[int]] = []
        for anchor in anchors:
            if not groups:
                groups.append([anchor]); continue
            previous = groups[-1][-1]
            between = signals[previous + 1:anchor]
            incompatible = any(signal and signal != module for signal in between)
            if not incompatible:
                groups[-1].append(anchor)
            else:
                groups.append([anchor])
        for group in groups:
            if len(group) < 2:
                continue
            start_i, end_i = group[0], group[-1]
            segment_proteins = ordered[start_i:end_i + 1]
            segment_id = f"module_{module}_{segment_proteins[0].start}_{segment_proteins[-1].end}"
            segments.append({"module_id": segment_id, "functional_module": module, "start": segment_proteins[0].start, "end": segment_proteins[-1].end, "protein_ids": [p.protein_id for p in segment_proteins], "informative_anchor_proteins": [ordered[i].protein_id for i in group], "contextual_unknown_proteins": [p.protein_id for i, p in enumerate(segment_proteins, start_i) if i not in group], "confidence": "HIGH", "boundary_status": "EXPLICIT_ANCHORS", "ambiguous_boundaries": False, "reasoning_summary": f"Two or more coherent {module} anchors are connected by compatible intervening proteins.", "rules_version": MODULE_RULES_VERSION, "_start_index": start_i, "_end_index": end_i})
    # Prevent overlapping segments from different labels; retain ambiguity on overlap.
    occupied: dict[int, dict[str, Any]] = {}
    canonical: list[dict[str, Any]] = []
    for segment in sorted(segments, key=lambda item: (item["_start_index"], item["_end_index"], item["functional_module"])):
        overlap = any(i in occupied for i in range(segment["_start_index"], segment["_end_index"] + 1))
        if overlap:
            segment["boundary_status"] = "AMBIGUOUS_OVERLAP"
            segment["ambiguous_boundaries"] = True
            continue
        canonical.append(segment)
        for i in range(segment["_start_index"], segment["_end_index"] + 1): occupied[i] = segment
    records: list[dict[str, Any]] = []
    for index, protein in enumerate(ordered):
        neighbors = ordered[max(0, index - window):index] + ordered[index + 1:index + window + 1]
        upstream = ordered[index - 1] if index else None
        downstream = ordered[index + 1] if index + 1 < len(ordered) else None
        def item(p):
            c = results.get(p.protein_id, {})
            return {"protein_id": p.protein_id, "start": p.start, "end": p.end, "strand": p.strand, "functional_state": c.get("functional_state"), "proposed_function": c.get("proposed_function"), "functional_category": c.get("functional_category"), "conservation_status": c.get("conservation_status")}
        neighbor_items = [item(p) for p in neighbors]
        containing = [segment for segment in canonical if segment["_start_index"] <= index <= segment["_end_index"]]
        module_record = containing[0] if containing else None
        local_signals = [s for s in signals[max(0, index-window):index+window+1] if s]
        ambiguity = len(set(local_signals)) > 1 and module_record is None
        module = module_record["functional_module"] if module_record else None
        module_start = module_record["start"] if module_record else None
        module_end = module_record["end"] if module_record else None
        records.append({"protein_id": protein.protein_id, "gene_order_index": index, "start": protein.start, "end": protein.end, "strand": protein.strand, "upstream_protein_id": upstream.protein_id if upstream else None, "downstream_protein_id": downstream.protein_id if downstream else None, "upstream_intergenic_distance": protein.start - upstream.end - 1 if upstream else None, "downstream_intergenic_distance": downstream.start - protein.end - 1 if downstream else None, "neighborhood_window": window, "neighbors": neighbor_items, "neighboring_proposed_functions": [n["proposed_function"] for n in neighbor_items if n["proposed_function"]], "neighboring_functional_categories": [n["functional_category"] for n in neighbor_items if n["functional_category"]], "neighboring_functional_states": [n["functional_state"] for n in neighbor_items if n["functional_state"]], "neighboring_conservation_states": [n["conservation_status"] for n in neighbor_items if n["conservation_status"]], "module_id": module_record["module_id"] if module_record else None, "functional_module": module, "module_start": module_start, "module_end": module_end, "context_confidence": "HIGH" if module_record else "LOW" if ambiguity else "NONE", "ambiguous_module_boundary": ambiguity or bool(module_record and module_record["ambiguous_boundaries"]), "module_boundary_reason": "mixed neighboring module signals" if ambiguity else None, "functional_state": results.get(protein.protein_id, {}).get("functional_state"), "proposed_function": results.get(protein.protein_id, {}).get("proposed_function"), "functional_category": results.get(protein.protein_id, {}).get("functional_category"), "conservation_status": results.get(protein.protein_id, {}).get("conservation_status"), "conserved_unknown": results.get(protein.protein_id, {}).get("functional_state") == "CONSERVED_UNKNOWN", "rules_version": MODULE_RULES_VERSION})
    for segment in canonical:
        segment.pop("_start_index", None); segment.pop("_end_index", None)
    return records, canonical


def write_context(root: str | Path, records: list[dict[str, Any]], modules: list[dict[str, Any]]) -> None:
    root = Path(root)
    (root / "genomic_context.json").write_text(json.dumps(records, indent=2, sort_keys=True))
    (root / "modules.json").write_text(json.dumps(modules, indent=2, sort_keys=True))
    columns = ["protein_id", "gene_order_index", "start", "end", "strand", "upstream_protein_id", "downstream_protein_id", "upstream_intergenic_distance", "downstream_intergenic_distance", "module_id", "functional_module", "module_start", "module_end", "functional_state", "proposed_function", "functional_category", "conservation_status", "conserved_unknown", "ambiguous_module_boundary", "context_confidence"]
    with (root / "genomic_context.tsv").open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=columns, delimiter="\t"); w.writeheader(); w.writerows({k: r.get(k) for k in columns} for r in records)
    mcols = ["module_id", "functional_module", "start", "end", "protein_ids", "confidence", "ambiguous_boundaries", "rules_version"]
    with (root / "modules.tsv").open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=mcols, delimiter="\t"); w.writeheader(); w.writerows({**{k: m.get(k) for k in mcols}, "protein_ids": ";".join(m["protein_ids"])} for m in modules)
