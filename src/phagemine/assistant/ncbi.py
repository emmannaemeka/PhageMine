"""NCBI submission support: error triage and guarded correction previews.

PhageMine never uploads anything to NCBI. This module only:

* classifies pasted validator / GenBank-staff messages into three outcomes;
* validates proposed corrections against what ``phagemine revise`` may change;
* refuses product names more specific than the recorded evidence supports;
* refuses metadata values the researcher never typed.

The classifier is keyword-based. Extend ``CATALOGUE`` with real reports from
past submissions; unknown messages are reported verbatim, never guessed at.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from phagemine.models import SubmissionMetadata

REVISABLE_FIELDS = {"product", "note", "partial", "locus_tag"}
METADATA_FIELDS = {f for f in SubmissionMetadata.__dataclass_fields__}

NEEDS_REANALYSIS = "NEEDS_REANALYSIS"
FIXABLE_LABEL = "FIXABLE_LABEL"
FIXABLE_METADATA = "FIXABLE_METADATA"
UNKNOWN = "UNKNOWN"


@dataclass(frozen=True)
class Rule:
    outcome: str
    pattern: str
    guidance: str


CATALOGUE: tuple[Rule, ...] = (
    Rule(NEEDS_REANALYSIS, r"internal\s*stop|internalstop", "A CDS translation contains a stop codon. The gene model or sequence must be checked; labels cannot fix this."),
    Rule(NEEDS_REANALYSIS, r"frame\s*-?\s*shift|frameshift", "Possible frameshift. Check the assembly at this locus and re-run annotation."),
    Rule(NEEDS_REANALYSIS, r"\bno\s*stop\b|nostop|missing stop|stop codon", "Stop codon problem in a gene model. Requires gene-model review or re-analysis."),
    Rule(NEEDS_REANALYSIS, r"start\s*codon|startcodon", "Start codon problem. Requires gene-model review."),
    Rule(NEEDS_REANALYSIS, r"translen|translation length|mismatchaa|translation.*(differ|mismatch)", "Translation does not match the CDS. Requires gene-model review."),
    Rule(NEEDS_REANALYSIS, r"overlap", "Overlapping features. Check whether both gene models are real; coordinates cannot be changed by revision."),
    Rule(NEEDS_REANALYSIS, r"short\s*(protein|cds|exon|prot)|shortprot", "Very short CDS. Review whether the gene model should be kept."),
    Rule(NEEDS_REANALYSIS, r"\bgap\b|beginsorendsingap", "Feature touches a sequence gap. Check the assembly."),
    Rule(FIXABLE_METADATA, r"collection[\s_-]*date|country|geo[\s_-]*loc|\bhost\b|isolation[\s_-]*source|biosource|lineage|\borganism\b|source\s*modifier|\bstrain\b|\bisolate\b|bioproject|biosample|submitter|author", "Submission metadata. Supply the correct value; PhageMine will not invent it."),
    Rule(FIXABLE_LABEL, r"product\s*name|productname|suspect\s*product|hypothetical|\bproduct\b", "Product name wording. Can be revised, but never made more specific than the evidence."),
    Rule(FIXABLE_LABEL, r"locus[\s_-]*tag", "Locus tag. Can be revised."),
    Rule(FIXABLE_LABEL, r"\bpartial", "Partial flag. Revisable, but check the gene ends first: a wrong partial flag can hide a coordinate problem."),
    Rule(FIXABLE_LABEL, r"\bnote\b", "Feature note. Can be revised."),
)

_SIGNAL = re.compile(r"error|warning|reject|fatal|problem|please|invalid|missing|suspect|\bseq_\w+", re.I)
_PROTEIN = re.compile(r"\bPM_\d{6}\b")
_LOCUS = re.compile(r"\b[A-Z][A-Z0-9]{2,11}_\d{3,6}\b")


def classify(text: str, limit: int = 60) -> dict:
    """Split pasted text into messages and assign each an outcome."""
    lines = [line.strip() for line in (text or "").splitlines() if line.strip()]
    items = []
    for line in lines:
        rule = next((r for r in CATALOGUE if re.search(r.pattern, line, re.I)), None)
        if rule is None and not _SIGNAL.search(line):
            continue
        items.append({"message": line[:400],
                      "outcome": rule.outcome if rule else UNKNOWN,
                      "guidance": rule.guidance if rule else "Not in PhageMine's error catalogue. Shown as received; no fix proposed.",
                      "protein_ids": sorted(set(_PROTEIN.findall(line))),
                      "locus_tags": sorted(set(_LOCUS.findall(line)) - set(_PROTEIN.findall(line)))})
    counts: dict[str, int] = {}
    for item in items:
        counts[item["outcome"]] = counts.get(item["outcome"], 0) + 1
    return {"messages": items[:limit], "total_messages": len(items), "counts": counts,
            "method": "keyword catalogue; unknown messages are not interpreted"}


def _norm(text: str | None) -> str:
    return re.sub(r"[^a-z0-9]+", " ", (text or "").lower()).strip()


def evidence_texts(evidence_row: dict, annotation_row: dict) -> list[str]:
    texts = [annotation_row.get("product"), annotation_row.get("proposed_function"),
             evidence_row.get("annotation")]
    for item in evidence_row.get("evidence", []) or []:
        if item.get("supports"):
            texts += [item.get("statement"), item.get("description"), item.get("family_name")]
    return [t for t in texts if t]


def check_product(new: str, current: str | None, supporting: list[str]) -> str | None:
    """Return a reason to refuse, or None when the new product name is allowed."""
    candidate = _norm(new)
    if not candidate:
        return "product name is empty"
    if candidate == "hypothetical protein":
        return None
    if candidate == _norm(current):
        return None
    if current and len(candidate) >= 4 and candidate in _norm(current):
        return None  # a less specific form of the current name
    if any(candidate == _norm(t) or (len(candidate) >= 6 and candidate in _norm(t)) for t in supporting):
        return None
    return ("more specific than the recorded evidence supports; allowed choices are "
            "'hypothetical protein', the current name, or wording found in this protein's evidence")


def check_metadata(metadata: dict, user_text: str) -> list[str]:
    problems = []
    corpus = user_text.lower()
    for key, value in metadata.items():
        if key not in METADATA_FIELDS:
            problems.append(f"{key} is not a submission metadata field")
            continue
        values = value if isinstance(value, list) else [value]
        for item in values:
            if not isinstance(item, str) or not item.strip():
                problems.append(f"{key} must be text")
            elif item.strip().lower() not in corpus:
                problems.append(f"{key} value '{item}' was not provided by the researcher in this conversation; ask for it")
    return problems
