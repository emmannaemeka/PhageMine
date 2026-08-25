"""Deterministic, conservative fusion of acquired evidence."""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path
from typing import Any

from .models import Evidence, Protein

FUSION_RULES_VERSION = "1.4"
EVIDENCE_HIERARCHY_VERSION = "1.1"
DIAGNOSTIC_DOMAIN_RULES_VERSION = "1.0"
DISPLAY_STATES = {
    "KNOWN_FUNCTION": "Specific function strongly supported",
    "PROBABLE_FUNCTION": "Likely function supported by evidence",
    "FUNCTIONAL_CLASS_ONLY": "Protein domain detected; full function unknown",
    "CONSERVED_UNKNOWN": "Conserved in phages; function unknown",
    "UNRESOLVED": "No reliable function identified",
    "CONFLICTING_EVIDENCE": "Conflicting evidence; review required",
}
MODALITY_FAMILIES = {
    "Pfam": "DOMAIN", "domain": "DOMAIN",
    "Swiss-Prot": "CURATED_SEQUENCE_HOMOLOGY", "sequence_homology": "CURATED_SEQUENCE_HOMOLOGY",
    "VOGDB": "VIRAL_ORTHOLOGY", "viral_orthology": "VIRAL_ORTHOLOGY",
    "PHROGs": "PHAGE_ORTHOLOGY", "phage_orthology": "PHAGE_ORTHOLOGY",
}
INCOMPATIBLE_FUNCTION_GROUPS = (
    {"integrase", "recombinase", "excisionase"},
    {"capsid", "portal protein", "terminase", "tail protein"},
    {"holin", "endolysin", "lysis"},
)
UNKNOWN_LABELS = {
    "", "hypothetical", "hypothetical protein", "unknown protein", "unknown function",
    "protein of unknown function", "uncharacterized protein", "uncharacterised protein",
    "conserved hypothetical protein", "conserved protein of unknown function",
    "putative uncharacterized protein", "no annotation", "none", "null",
}

# These descriptions can be legitimate database similarities, but they are not
# defensible bacteriophage product names without independent, phage-specific
# corroboration.  They remain visible in the raw evidence output.
TAXON_SPECIFIC_TERMS = (
    "interferon", "apolipoprotein", "kinetochore", "centromere", "exocyst",
    "male sterility", "cerebellar degeneration", "baculovirus", "mitochondrial",
    "chloroplast", "chloroplastic", "arabidopsis", "human ", "mouse ",
    "anthrax toxin", "lethal factor", "ino80", "balf1", "arabinogalactan",
    "fungal", "apoptosis", "outer membrane lipoprotein",
    "spore coat", "starch initiation", "znf598", "plant ", "animal ",
    "queuosine salvage protein",
)
PHAGE_SAFE_FUNCTION_TERMS = (
    "phage", "virion", "capsid", "portal", "terminase", "tail", "baseplate",
    "holin", "endolysin", "lysin", "integrase", "recombinase", "excisionase",
    "polymerase", "helicase", "nuclease", "methyltransferase", "atpase",
    "peptidase", "glycosyl hydrolase", "dna-binding", "rna-binding",
    "thymidylate synthase", "ribonucleotide reductase", "scaffolding protein",
)


def _independently_corroborated(label: str, evidence: Evidence, accepted: list[Evidence]) -> bool:
    for other in accepted:
        if other is evidence or other.source == evidence.source or other.source == "Pfam":
            continue
        other_label = normalize_function(other.description)
        if other_label and (label == other_label or label in other_label or other_label in label):
            return True
    return False


def _product_candidate(evidence: Evidence, label: str | None, accepted: list[Evidence]) -> tuple[str | None, str | None]:
    """Return a defensible whole-protein product candidate.

    Pfam establishes domain architecture, not a complete protein product.
    Organism-specific descriptions remain evidence but cannot be transferred to
    a phage product field.  Existing conservative family-level rewrites are
    applied after these source-level rules.
    """
    if not label:
        return None, None
    if evidence.source == "Pfam" or evidence.modality == "domain":
        return None, "domain-only evidence retained as a note, not transferred as a protein product"
    conservative, flag = _conservative_label(label, accepted)
    # A supported family-level rewrite (for example AAA-family ATPase) is
    # useful even when the matched database member had an unsafe organelle-
    # specific name.
    if conservative and conservative != label:
        return conservative, flag
    if any(term in label for term in TAXON_SPECIFIC_TERMS):
        return None, "taxon-specific product label suppressed because it lacks independent phage-specific support"
    safe_function = any(term in label for term in PHAGE_SAFE_FUNCTION_TERMS)
    corroborated = _independently_corroborated(label, evidence, accepted)
    if evidence.source == "VOGDB" and not (safe_function or corroborated):
        return None, "VOG member description retained as evidence but not transferred without a phage-compatible function or independent corroboration"
    if evidence.source == "Swiss-Prot":
        organism = str(evidence.metrics.get("organism") or "").lower()
        phage_record = "phage" in organism or "bacteriophage" in organism
        if not (phage_record or safe_function or corroborated):
            return None, "Swiss-Prot product retained as homology evidence but not transferred without phage-compatible provenance or corroboration"
    return conservative, flag


def _domain_summary(accepted: list[Evidence]) -> str | None:
    labels = []
    for evidence in accepted:
        if evidence.source != "Pfam" and evidence.modality != "domain":
            continue
        label = normalize_function(evidence.description)
        if label and "phage ninh protein and transposase" in label:
            label = "ninh-like domain"
        elif label and "enterobacter phage enc34, ssdna-binding protein" in label:
            label = "single-stranded dna-binding domain"
        elif label and "disarm protein drme" in label:
            label = "atpase-related domain"
        if label and not any(term in label for term in TAXON_SPECIFIC_TERMS) and label not in labels:
            labels.append(label)
    if any("ninh" in label for label in labels):
        labels = [label for label in labels if "transposase" not in label and "ninh" not in label]
        labels.insert(0, "ninh-like domain")
    if any("atpase-related domain" in label for label in labels):
        labels = [label for label in labels if "p-type atpase actuator" not in label and "atpase-related domain" not in label]
        labels.insert(0, "atpase-related domain")
    if not labels:
        return None
    shown = labels[:2]
    suffix = " and ".join(shown)
    return f"hypothetical protein containing {suffix}"


def normalize_function(description: str | None) -> str | None:
    if description is None:
        return None
    value = re.sub(r"\s+", " ", description.strip().lower())
    value = re.sub(r"^refseq\s+", "", value)
    value = re.sub(r"^sp\|[^|]+\|", "", value)
    # Swiss-Prot descriptions sometimes retain an entry name (for example
    # ``YR614_MIMIV``) after the accession has been removed.  Entry/locus names
    # identify a database record, not a transferable biological function.
    value = re.sub(r"^[a-z0-9]+_[a-z0-9]+\s+", "", value)
    value = re.sub(r"\s*\{eco:[^}]+\}\s*$", "", value, flags=re.IGNORECASE).strip()
    canonical = {
        "major head protein": "major capsid protein",
        "hoc-like head decoration": "hoc-like head decoration protein",
        "head scaffolding protein": "head scaffolding protein",
        "head closure hc1": "head closure protein hc1",
        "dna helicase": "dna helicase",
        "terminase large subunit": "terminase large subunit",
        "endolysin": "endolysin",
        "tail spike": "tail spike protein",
        "head morphogenesis": "head morphogenesis protein",
        "terminase": "terminase protein",
    }
    value = canonical.get(value, value)
    return None if value in UNKNOWN_LABELS else value


def _numeric_metric(evidence: Evidence, *names: str) -> float | None:
    for name in names:
        value = evidence.metrics.get(name)
        try:
            if value not in (None, "", "No_PHROG", "No_PHROGs_HMM"):
                return float(value)
        except (TypeError, ValueError):
            pass
    return None


def _candidate_score(evidence: Evidence, label: str) -> float:
    """Rank product hypotheses using provenance and alignment support."""
    source = {"Swiss-Prot": 50.0, "PHROGs": 45.0, "VOGDB": 30.0}.get(evidence.source, 10.0)
    strength = {"EXPERIMENTAL": 40.0, "STRONG": 25.0, "MODERATE": 12.0, "WEAK": 0.0}.get(evidence.evidence_strength, 0.0)
    identity = _numeric_metric(evidence, "percent_identity", "sequence_identity", "identity")
    if identity is not None and identity > 1:
        identity /= 100.0
    qcov = _numeric_metric(evidence, "query_coverage", "qcov")
    if qcov is not None and qcov > 1:
        qcov /= 100.0
    bits = _numeric_metric(evidence, "bit_score", "bitscore", "mmseqs_score", "score")
    evalue = _numeric_metric(evidence, "evalue")
    score = source + strength
    score += 20.0 * max(0.0, min(identity or 0.0, 1.0))
    score += 15.0 * max(0.0, min(qcov or 0.0, 1.0))
    score += min(max(bits or 0.0, 0.0) / 50.0, 12.0)
    if evalue is not None and evalue <= 1e-20:
        score += 8.0
    if any(term in label for term in PHAGE_SAFE_FUNCTION_TERMS):
        score += 8.0
    return score


def _select_product(informative: list[tuple[int, Evidence, str]]) -> tuple[str | None, Evidence | None, list[dict[str, Any]]]:
    """Select one auditable product or abstain when evidence is unresolved."""
    ranked = sorted(
        ((i, evidence, label, _candidate_score(evidence, label)) for i, evidence, label in informative),
        key=lambda item: (-item[3], item[2], item[0]),
    )
    alternatives = [
        {"label": label, "source": evidence.source,
         "identifier": evidence.identifier or evidence.family_name,
         "adjudication_score": round(score, 3)}
        for _, evidence, label, score in ranked
    ]
    if not ranked:
        return None, None, alternatives
    _, winner, label, winner_score = ranked[0]
    # Prefer a supported specific phage role over its own broad parent term.
    # This is a terminology hierarchy, not a score override between unrelated
    # biological functions.
    specific = {
        "head protein": "head morphogenesis protein",
        "nuclease": "holliday junction resolvase",
        "endonuclease": "holliday junction resolvase",
    }
    labels_present = {item[2] for item in ranked}
    preferred = specific.get(label)
    if preferred and preferred in labels_present:
        chosen = next(item for item in ranked if item[2] == preferred)
        _, winner, label, winner_score = chosen
        return label, winner, alternatives
    runner_score = ranked[1][3] if len(ranked) > 1 else float("-inf")
    overlapping = all(label == item[2] or label in item[2] or item[2] in label for item in ranked[1:])
    quantified = any(_numeric_metric(winner, name) is not None for name in (
        "percent_identity", "sequence_identity", "identity", "query_coverage",
        "qcov", "bit_score", "bitscore", "mmseqs_score", "score", "evalue",
    ))
    if len(ranked) == 1 or overlapping or (quantified and winner_score - runner_score >= 8.0):
        return label, winner, alternatives
    return None, None, alternatives


def _diagnostic_product(accepted: list[Evidence]) -> tuple[str | None, Evidence | None, str | None]:
    """Promote only narrowly reviewed domain/orthology combinations.

    This is deliberately a small allow-list.  It prevents arbitrary Pfam text
    from becoming a product while allowing diagnostic enzyme families and
    compatible phage orthology to yield useful, qualified annotations.
    """
    records = [(e, normalize_function(e.description)) for e in accepted]
    records = [(e, label) for e, label in records if label]
    pfam = [(e, label) for e, label in records if e.source == "Pfam" or e.modality == "domain"]
    phrogs = [(e, label) for e, label in records if e.source == "PHROGs"]
    categories = " ".join(filter(None, (_category(e) for e in accepted))).lower()
    domain_text = " ".join(label for _, label in pfam)

    if "dna polymerase family a" in domain_text:
        evidence = next(e for e, label in pfam if "dna polymerase family a" in label)
        return "family-a dna polymerase", evidence, "diagnostic DNA polymerase family A domain"
    if "enterobacter phage enc34, ssdna-binding protein" in domain_text or "single-stranded dna-binding" in domain_text:
        evidence = next(e for e, label in pfam if "ssdna-binding" in label or "single-stranded dna-binding" in label)
        return "single-stranded dna-binding protein", evidence, "diagnostic phage ssDNA-binding domain"
    if "vrr-nuc" in domain_text:
        resolvase = next(((e, label) for e, label in phrogs if "holliday junction resolvase" in label), None)
        if resolvase:
            return "holliday junction resolvase", resolvase[0], "VRR-Nuc domain corroborated by PHROGs resolvase orthology"
    hth = "helix-turn-helix" in domain_text
    excisionase = next(((e, label) for e, label in phrogs if "excisionase" in label and "transcriptional regulator" in label), None)
    if hth and excisionase and "integration and excision" in categories:
        return "excisionase and transcriptional regulator", excisionase[0], "HTH domains and integration/excision category corroborate PHROGs orthology"
    regulator = next(((e, label) for e, label in records if "transcriptional regulator" in label), None)
    if hth and regulator:
        return "helix-turn-helix transcriptional regulator", regulator[0], "HTH domains corroborate a broad transcriptional-regulator assignment"
    return None, None, None


def _gene_and_ec(accepted: list[Evidence]) -> tuple[str | None, str | None]:
    """Transfer identifiers only from strong, accepted, explicit records."""
    genes: list[str] = []; ecs: list[str] = []
    for evidence in accepted:
        if evidence.evidence_strength not in {"STRONG", "EXPERIMENTAL"}:
            continue
        gene = evidence.metrics.get("gene") or evidence.metrics.get("gene_name")
        ec = evidence.metrics.get("ec_number") or evidence.metrics.get("ec")
        if gene and re.fullmatch(r"[A-Za-z][A-Za-z0-9_.-]{0,19}", str(gene)):
            genes.append(str(gene))
        if ec:
            match = re.search(r"\b\d+\.(?:\d+|-)\.(?:\d+|-)\.(?:\d+|-)\b", str(ec))
            if match:
                ecs.append(match.group(0))
    return (genes[0] if genes and len(set(genes)) == 1 else None,
            ecs[0] if ecs and len(set(ecs)) == 1 else None)


def _conservative_label(label: str | None, accepted: list[Evidence]) -> tuple[str | None, str | None]:
    """Prevent taxon-, organelle-, and locus-specific over-annotation.

    A significant homology/orthology match establishes relatedness.  It does
    not, by itself, establish that a phage protein performs the exact role of a
    named mitochondrial, chloroplast, or eukaryotic-virus database member.
    Return a defensible family-level label and an audit flag when independent
    domain evidence permits one; otherwise suppress the unsafe specific label.
    """
    if not label:
        return None, None
    evidence_text = " ".join(
        filter(None, (normalize_function(item.description) for item in accepted))
    )
    if "bcs1" in label and any(term in evidence_text for term in ("aaa", "atpase")):
        return "bcs1-like aaa-family atpase", "organelle-specific BCS1 label reduced to a family-level ATPase annotation"
    if "band 7" in label or "spfh" in label:
        return "band 7/spfh family protein", "record-specific Band 7 label reduced to a family-level annotation"
    if re.search(r"\bopg\d+\b", label) or "immune evasion protein opg" in label:
        if "kelch" in evidence_text or "beta-propeller" in evidence_text:
            return "kelch-repeat beta-propeller protein", "eukaryotic-virus locus/function label reduced to the independently supported domain architecture"
        return None, "eukaryotic-virus locus/function label suppressed because its specific function lacks independent support"
    if any(term in label for term in ("mitochondrial", "chloroplastic", "chloroplast")):
        return None, "organelle-specific product label suppressed because its biological context lacks independent support"
    return label, None


def _category(evidence: Evidence) -> str | None:
    value = evidence.metrics.get("functional_category") or evidence.metrics.get("category")
    if value is None or str(value).strip().lower() in UNKNOWN_LABELS | {"xu", "unknown"}:
        return None
    return str(value).strip()


def _conservation(evidence: list[Evidence]) -> str:
    accepted = [e for e in evidence if e.supports and e.source in {"PHROGs", "VOGDB"}]
    strong = [e for e in accepted if e.evidence_strength == "STRONG"]
    if strong:
        return "STRONGLY_CONSERVED"
    if accepted:
        return "CONSERVED"
    return "NOT_ESTABLISHED"


def classify_protein(protein: Protein) -> dict[str, Any]:
    accepted = [e for e in protein.evidence if e.supports]
    gene, ec_number = _gene_and_ec(accepted)
    conservative_flags: list[str] = []
    conservative_rewrites: list[str] = []
    informative = []
    for i, evidence in enumerate(accepted):
        original_label = normalize_function(evidence.description)
        label, flag = _product_candidate(evidence, original_label, accepted)
        if flag and flag not in conservative_flags:
            conservative_flags.append(flag)
        if flag and label and label != original_label:
            conservative_rewrites.append(label)
        if label:
            informative.append((i, evidence, label))
    strong_info = [(i, e, label) for i, e, label in informative if e.evidence_strength == "STRONG"]
    categories = sorted({_category(e) for e in accepted if _category(e)})
    support_ids = [f"{e.source}:{e.identifier or e.family_name or i}" for i, e in enumerate(accepted)]
    ortholog_groups = sorted({
        str(e.identifier or e.family_name) for e in accepted
        if e.source in {"PHROGs", "VOGDB"} and (e.identifier or e.family_name)
    })
    supporting_sources = sorted({e.source for e in accepted})
    supporting_modalities = sorted({MODALITY_FAMILIES.get(e.source, MODALITY_FAMILIES.get(e.modality, e.modality)) for e in accepted})
    conflicting = []
    labels_by_source: dict[str, set[str]] = {}
    for _, e, label in strong_info:
        labels_by_source.setdefault(e.source, set()).add(label)
    strong_labels = {label for _, _, label in strong_info}
    compatible = all(a == b or a in b or b in a for a in strong_labels for b in strong_labels)
    explicit_incompatibility = any(any(any(term in label for term in group) for label in strong_labels) and
                                   any(any(term in label for term in other) for label in strong_labels)
                                   for index, group in enumerate(INCOMPATIBLE_FUNCTION_GROUPS)
                                   for other in INCOMPATIBLE_FUNCTION_GROUPS[index + 1:])
    ambiguity = list(conservative_flags)
    if len(strong_labels) > 1 and not compatible and not explicit_incompatibility:
        ambiguity.append("distinct strong labels could not be confidently established as biologically incompatible")
    if len(strong_labels) > 1 and len(labels_by_source) > 1 and explicit_incompatibility:
        conflicting = [(i, e, label) for i, e, label in strong_info]
    conflict_ids = [f"{e.source}:{e.identifier or e.family_name or i}" for i, e, _ in conflicting]
    conflict_sources = sorted({e.source for _, e, _ in conflicting})
    conflict_descriptions = [e.description for _, e, _ in conflicting if e.description]
    labels = sorted({label for _, _, label in informative})
    rewrite_labels = sorted(set(conservative_rewrites))
    # A safety rewrite is deliberately more conservative than the matched
    # member name and therefore takes precedence over ancillary domain labels.
    selected_product, selected_evidence, product_alternatives = _select_product(informative)
    diagnostic_product, diagnostic_evidence, diagnostic_reason = _diagnostic_product(accepted)
    conflict_resolved = bool(conflicting and diagnostic_product)
    if diagnostic_product:
        selected_product, selected_evidence = diagnostic_product, diagnostic_evidence
        if diagnostic_reason and diagnostic_reason not in ambiguity:
            ambiguity.append(diagnostic_reason)
    elif conflicting:
        selected_product, selected_evidence = None, None
    if len(rewrite_labels) == 1:
        selected_product = rewrite_labels[0]
        selected_evidence = next((e for _, e, label in informative if label == selected_product), selected_evidence)
    domain_summary = _domain_summary(accepted)
    conservation = _conservation(protein.evidence)
    display_product = selected_product or domain_summary or ("conserved phage protein of unknown function" if conservation in {"STRONGLY_CONSERVED", "CONSERVED"} else "hypothetical protein")
    independent_strong = len({e.source for _, e, _ in strong_info}) >= 2
    swiss_strong = any(e.source == "Swiss-Prot" for _, e, _ in strong_info)
    strong_families = {MODALITY_FAMILIES.get(e.source, MODALITY_FAMILIES.get(e.modality, e.modality)) for _, e, _ in strong_info}
    orthology_only = strong_families and strong_families <= {"PHAGE_ORTHOLOGY", "VIRAL_ORTHOLOGY"}
    agreeing_sources = len({e.source for _, e, label in informative if selected_product and label == selected_product})
    if conflicting and not conflict_resolved:
        state, confidence, reasons = "CONFLICTING_EVIDENCE", "LOW", ["accepted STRONG informative evidence supports incompatible normalized functions"]
    elif strong_info and (independent_strong or (swiss_strong and agreeing_sources >= 2)) and len(strong_families) >= 2 and not orthology_only and not ambiguity:
        state, confidence = "KNOWN_FUNCTION", "HIGH"
        reasons = [f"STRONG {e.source} informative evidence" for _, e, _ in strong_info]
        if independent_strong:
            reasons.append("independent agreeing source corroboration")
        reasons.append("no accepted strong conflict")
    elif selected_product:
        low_diagnostic = diagnostic_reason in {
            "diagnostic phage ssDNA-binding domain",
            "HTH domains corroborate a broad transcriptional-regulator assignment",
        }
        state = "PROBABLE_FUNCTION"
        confidence = "LOW" if low_diagnostic else ("MODERATE" if strong_info or diagnostic_product else "LOW")
        reasons = [
            f"accepted {e.evidence_strength or 'informative'} {e.source} evidence supports a functional interpretation"
            for _, e, _ in informative
        ]
        if diagnostic_reason:
            reasons.append(diagnostic_reason)
    elif categories or domain_summary or any(e.source == "Pfam" or e.modality == "domain" for e in accepted):
        state, confidence = "FUNCTIONAL_CLASS_ONLY", "LOW"
        reasons = ["accepted evidence supports domain architecture or a broad functional category, but not a complete protein function"]
    elif conservation in {"STRONGLY_CONSERVED", "CONSERVED"}:
        state, confidence = "CONSERVED_UNKNOWN", "LOW"
        reasons = ["accepted PHROGs/VOGDB orthology supports conservation but no accepted informative functional evidence is present"]
    else:
        state, confidence, reasons = "UNRESOLVED", "NONE", ["no accepted evidence establishes a function, category, or meaningful conservation"]
    if selected_product:
        interpretation = "Defensible computational product hypothesis; review supporting alignments before biological assertion."
    elif domain_summary:
        interpretation = "A conserved domain was detected, but it does not establish the complete protein function."
    elif conservation in {"STRONGLY_CONSERVED", "CONSERVED"}:
        interpretation = "Viral/phage orthology supports conservation, but the biological function remains unknown."
    else:
        interpretation = "Predicted coding sequence with no accepted functional or conservation evidence; function remains unknown."
    strength_rank = {"EXPERIMENTAL": 4, "STRONG": 3, "MODERATE": 2, "WEAK": 1, None: 0}
    source_rank = {"PHROGs": 4, "VOGDB": 3, "Swiss-Prot": 2, "Pfam": 1}
    ranked_evidence = sorted(
        accepted,
        key=lambda e: (-strength_rank.get(e.evidence_strength, 0), -source_rank.get(e.source, 0), str(e.identifier or e.family_name or "")),
    )
    matching_phrogs = next(
        (e for _, e, label in informative if selected_product and label == selected_product and e.source == "PHROGs"),
        None,
    )
    best = matching_phrogs or selected_evidence or (ranked_evidence[0] if ranked_evidence else None)
    best_evidence = (
        f"{best.source}:{best.identifier or best.family_name or 'match'}"
        + (f" — {best.description}" if best.description else "")
        if best else "No accepted evidence"
    )
    enzyme_domain = any(term in " ".join(filter(None, (normalize_function(e.description) for e in accepted if e.source == "Pfam" or e.modality == "domain")))
                        for term in ("atpase", "polymerase", "helicase", "nuclease", "transposase"))
    if protein.length < 80 and enzyme_domain:
        review_flag = "POSSIBLE_PARTIAL_OR_FALSE_ORF"
    else:
        review_flag = "REVIEW_REQUIRED" if state in {"CONFLICTING_EVIDENCE", "UNRESOLVED"} else "NONE"
    experimental = any(e.evidence_strength == "EXPERIMENTAL" or e.level.value == "experimentally established" for e in accepted)
    if experimental and selected_product:
        evidence_tier, evidence_tier_label = 1, "Experimentally supported database evidence"
    elif state == "KNOWN_FUNCTION":
        evidence_tier, evidence_tier_label = 2, "Curated and independently corroborated function"
    elif state == "PROBABLE_FUNCTION" and len(supporting_sources) >= 2:
        evidence_tier, evidence_tier_label = 3, "Independent computational sources support the same function"
    elif state in {"PROBABLE_FUNCTION", "FUNCTIONAL_CLASS_ONLY"}:
        evidence_tier, evidence_tier_label = 4, "Single-source function or protein-domain support"
    elif state == "CONSERVED_UNKNOWN":
        evidence_tier, evidence_tier_label = 5, "Phage conservation established; function unknown"
    else:
        evidence_tier, evidence_tier_label = 6, "Function unresolved or conflicting"
    return {
        "protein_id": protein.protein_id, "functional_state": state,
        "display_classification": DISPLAY_STATES[state],
        "proposed_function": selected_product, "normalized_function": selected_product,
        "display_product": display_product, "domain_summary": domain_summary,
        "gene": gene, "ec_number": ec_number,
        "biotechnology_relevance": "Potential capsid-display candidate; experimental confirmation required." if "hoc-like head decoration protein" in display_product else None,
        "scientific_interpretation": interpretation,
        "best_evidence": best_evidence, "review_flag": review_flag,
        "product_alternatives": product_alternatives,
        "diagnostic_domain_rule": diagnostic_reason,
        "diagnostic_domain_rules_version": DIAGNOSTIC_DOMAIN_RULES_VERSION,
        "conflict_resolved_by_corroboration": conflict_resolved,
        "evidence_tier": evidence_tier, "evidence_tier_label": evidence_tier_label,
        "evidence_hierarchy_version": EVIDENCE_HIERARCHY_VERSION,
        "functional_category": categories[0] if len(categories) == 1 else ("; ".join(categories) if categories else None),
        "confidence": confidence, "confidence_reasons": reasons,
        "conservation_status": conservation,
        "supporting_sources": supporting_sources, "supporting_modalities": supporting_modalities,
        "supporting_source_count": len(supporting_sources), "supporting_modality_count": len(supporting_modalities), "supporting_record_count": len(informative),
        "supporting_evidence_ids": support_ids,
        "ortholog_groups": ortholog_groups,
        "conflicting_sources": conflict_sources, "conflicting_evidence_ids": conflict_ids,
        "conflict_descriptions": conflict_descriptions, "ambiguity_flags": ambiguity,
        "reasoning_summary": "; ".join(reasons) + ".", "fusion_rules_version": FUSION_RULES_VERSION,
    }


def classify_proteins(proteins: list[Protein]) -> list[dict[str, Any]]:
    results = [classify_protein(protein) for protein in proteins]
    by_id = {item["protein_id"]: item for item in results}
    ordered = sorted(proteins, key=lambda p: (p.start, p.end, p.protein_id))

    def tail_signal(item: dict[str, Any]) -> bool:
        text = " ".join(str(item.get(key) or "") for key in
                        ("proposed_function", "functional_category")).lower()
        return any(term in text for term in
                   ("tail", "baseplate", "fiber", "spike", "neck", "connector"))

    for index, protein in enumerate(ordered):
        item = by_id[protein.protein_id]
        if item.get("proposed_function"):
            continue
        neighbors = ordered[max(0, index - 2):index] + ordered[index + 1:index + 3]
        anchors = [by_id[p.protein_id] for p in neighbors if tail_signal(by_id[p.protein_id])]
        if len(anchors) < 2 or item.get("conservation_status") not in {"CONSERVED", "STRONGLY_CONSERVED"}:
            continue
        alternatives = item.get("product_alternatives") or []
        minor_tail = next((candidate for candidate in alternatives
                           if candidate.get("label") == "minor tail protein"), None)
        domain_text = str(item.get("domain_summary") or "").lower()
        if minor_tail:
            product = "putative minor tail protein"
            basis = "accepted minor-tail candidate supported by a conserved tail-module neighbourhood"
        elif any(term in domain_text for term in ("tail tube", "tail protein", "mbg domain")) or item.get("functional_state") == "CONSERVED_UNKNOWN":
            product = "putative tail-associated protein"
            basis = "phage conservation and two or more neighbouring tail-module anchors; exact role remains unresolved"
        else:
            continue
        item.update({
            "proposed_function": product, "normalized_function": product,
            "display_product": product, "functional_state": "PROBABLE_FUNCTION",
            "display_classification": DISPLAY_STATES["PROBABLE_FUNCTION"],
            "confidence": "LOW", "review_flag": "NONE",
            "context_support": {
                "module": "tail", "neighboring_proteins": [anchor["protein_id"] for anchor in anchors],
                "neighboring_functions": [anchor.get("display_product") for anchor in anchors],
                "effect": f"upgraded to {product}", "context_rules_version": "1.0",
                "basis": basis,
            },
            "reasoning_summary": item.get("reasoning_summary", "") + " " + basis + ".",
        })
    return results


def write_classification(root: str | Path, proteins: list[Protein], results: list[dict[str, Any]] | None = None) -> None:
    results = results if results is not None else classify_proteins(proteins)
    by_id = {p.protein_id: p for p in proteins}
    path = Path(root)
    (path / "functional_classification.json").write_text(json.dumps(results, indent=2, sort_keys=True))
    columns = ["protein_id", "start", "end", "strand", "length_aa", "gene", "ec_number", "functional_state", "display_classification", "proposed_function", "display_product", "confidence", "evidence_tier", "evidence_tier_label", "best_evidence", "review_flag", "domain_summary", "functional_category", "biotechnology_relevance", "conservation_status", "ortholog_groups", "supporting_sources", "supporting_source_count", "supporting_record_count", "conflict", "scientific_interpretation", "reasoning_summary"]
    with (path / "functional_classification.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for result in results:
            protein = by_id[result["protein_id"]]
            writer.writerow({**{key: result.get(key) for key in columns}, "start": protein.start, "end": protein.end, "strand": protein.strand, "length_aa": protein.length, "ortholog_groups": ";".join(result["ortholog_groups"]), "supporting_sources": ";".join(result["supporting_sources"]), "conflict": bool(result["conflicting_evidence_ids"])})
