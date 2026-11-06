"""Offline, evidence-grounded request routing. Never evaluates generated code."""
from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass(frozen=True)
class Request:
    action: str
    query: str = ""


def interpret(text: str) -> Request:
    text = text.strip()
    if not text or len(text) > 1000:
        return Request("help")
    lower = text.lower().rstrip(".?!")
    if re.search(r"\b(?:delete|remove|overwrite|execute|shell|install)\b", lower):
        return Request("help")
    if re.fullmatch(r"(?:please )?(?:annotate|analyse|analyze)(?: (?:this|my|the|a))?(?: (?:phage|genome|file|fasta|genomes|files|uploaded|sequence|sequences))+", lower):
        return Request("annotate")
    if re.fullmatch(r"(?:please )?(?:summari[sz]e|explain)(?: (?:my|these|the))? (?:results|annotations)", lower):
        return Request("summary")
    if re.search(r"\b(?:compare|shared|common)\b", lower):
        return Request("compare")
    if re.fullmatch(r"(?:check |show )?(?:my |the )?(?:installation|environment|database status|databases)", lower):
        return Request("environment")
    match = re.fullmatch(r"(?:please )?(?:can you |could you )?(?:find|search for|show|explain|tell me about)\s+(?:me\s+)?(.+)", lower)
    if match:
        query = re.sub(r"\b(?:the|a|an|gene|genes|protein|proteins|responsible for|in my results|in these results)\b", " ", match[1])
        query = " ".join(query.split())
        if query and not re.search(r"\b(?:and|then|except|not|without)\b", query):
            return Request("search", query)
    return Request("help")


SEARCH_FIELDS = ("protein_id", "gene", "product", "proposed_function", "annotation", "functional_category")


def find_records(records: list[dict], query: str) -> list[dict]:
    terms = re.findall(r"[\w.-]+", query.lower())
    if not terms:
        return []
    return [row for row in records if all(term in " ".join(str(row.get(k) or "") for k in SEARCH_FIELDS).lower() for term in terms)]


def summarize(records: list[dict]) -> str:
    if not records:
        return "No annotation records were found. Select a completed result folder or start an annotation."
    samples = len({r.get("sample") for r in records})
    unknown = sum(any(word in str(r.get("product") or r.get("proposed_function") or "").lower()
                      for word in ("hypothetical", "unknown")) for r in records)
    conflicts = sum(bool(r.get("conflict_review")) for r in records)
    return (f"The selected outputs contain {len(records)} annotation records across {samples} sample(s). "
            f"{unknown} records have products explicitly labelled hypothetical or unknown; "
            f"{conflicts} have attached conflict-review records. "
            "These counts describe saved outputs, not experimentally verified functions. "
            "Confidence labels express evidence strength, not probabilities. Missing matches do not prove gene absence.")
