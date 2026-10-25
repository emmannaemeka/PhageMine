"""NCBI submission support: error triage and guarded correction previews.

PhageMine never uploads anything to NCBI. This module only:

* classifies pasted validator / GenBank-staff messages into three outcomes;
* validates proposed corrections against what ``phagemine revise`` may change;
* refuses product names more specific than PhageMine's own adjudicated name,
  including dropping uncertainty qualifiers such as "putative";
* refuses metadata values the researcher never typed.

Every correction still needs the researcher's explicit approval (see
``approvals``); these checks only limit what may be proposed.

The classifier is keyword-based. Extend ``CATALOGUE`` with real reports from
past submissions; unknown messages are reported verbatim, never guessed at.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from phagemine.models import SubmissionMetadata

REVISABLE_FIELDS = {"product", "note", "locus_tag"}
# Partial flags describe whether a gene's ends are complete: a curation decision, not wording.
CURATION_FIELDS = {"partial", "start", "end", "strand", "codon_start", "sequence", "translation"}
QUALIFIERS = ("putative", "probable", "possible", "predicted", "potential", "uncharacterized",
              "uncharacterised", "like", "domain containing", "family protein", "homolog")
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
    Rule(NEEDS_REANALYSIS, r"\bpartial", "Partial flag: whether the gene's ends are complete is a curation decision. Review the gene model before changing it."),
    Rule(FIXABLE_LABEL, r"\bnote\b", "Feature note. Can be revised."),
)

_PROTEIN = re.compile(r"\bPM_\d{6}\b")
_LOCUS = re.compile(r"\b[A-Z][A-Z0-9]{2,11}_\d{3,6}\b")


def classify(text: str, limit: int = 60) -> dict:
    """Split pasted text into messages and assign each an outcome."""
    # Keep every line with words in it; only separators and blank lines are skipped.
    lines = [line.strip() for line in (text or "").splitlines() if re.search(r"[A-Za-z0-9]", line)]
    items = []
    for line in lines:
        rule = next((r for r in CATALOGUE if re.search(r.pattern, line, re.I)), None)
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


def _qualifiers(name: str) -> tuple[set[str], str]:
    """Split a normalised product name into its uncertainty qualifiers and its core."""
    text, found = f" {_norm(name)} ", set()
    for word in QUALIFIERS:
        if f" {word} " in text:
            found.add(word)
            text = text.replace(f" {word} ", " ")
    return found, " ".join(text.split())


def check_product(new: str, current: str | None, _unused_evidence: list[str] | None = None) -> str | None:
    """Return a reason to refuse, or None when the new product name is allowed.

    Allowed: "hypothetical protein"; the current name re-formatted; or the current
    name with an added uncertainty qualifier. Anything more specific, including
    dropping "putative", needs curation with evidence, not an NCBI wording fix.
    Words found in evidence descriptions never authorise a name.
    """
    candidate = _norm(new)
    if not candidate:
        return "product name is empty"
    if candidate == "hypothetical protein" or candidate == _norm(current):
        return None
    if current:
        new_q, new_core = _qualifiers(new)
        old_q, old_core = _qualifiers(current)
        if new_core == old_core and old_q <= new_q:
            return None
        if new_core == old_core:
            return "it drops an uncertainty qualifier (" + ", ".join(sorted(old_q - new_q)) + ") from the current name"
    return ("more specific than PhageMine's adjudicated name; allowed choices are 'hypothetical protein', "
            "the current name, or the current name with an added qualifier such as 'putative'. "
            "A more specific name needs evidence review and curation")


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
