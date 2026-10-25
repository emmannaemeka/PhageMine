"""Deterministic checks that an answer only states what tools returned.

Evidence is tool output only; what the researcher typed is never evidence.
Checked: protein and job IDs, numbers, completion claims (bound to the job or
sample named in the same sentence), and the unknown-means-novel error.
Not checked: the meaning of every sentence. Evidence tables stay the authority.
"""
from __future__ import annotations

import json
import re

from .store import Project

PROTEIN_ID = re.compile(r"\bPM_\d{6}\b")
JOB_ID = re.compile(r"\b(?:annotation|comparison|ncbirevision)-[0-9a-f]{8}\b")
COMPLETION = re.compile(r"\b(?:completed?|finished|succeeded|successful(?:ly)?|is done|are done|has been (?:annotated|compared))\b", re.I)
HEDGE = re.compile(r"\b(?:not|n't|will|once|when|after|until|yet|still|running|queued|waiting|failed|if)\b", re.I)
NOVEL = re.compile(r"\bnovel\b", re.I)
NEGATION = re.compile(r"(?:not|n't|no|never|does not mean|doesn't mean|without)\W+(?:\w+\W+){0,4}$", re.I)
NUMBER = re.compile(r"(?<![\w.])(\d[\d,]*(?:\.\d+)?(?:[eE]-?\d+)?)(?![\w])")
ISO_DATE = re.compile(r"\b\d{4}-\d{2}-\d{2}\b")


def _sentences(text: str) -> list[str]:
    return [s for s in re.split(r"(?<=[.!?\n])\s+", text) if s.strip()]


def _numbers(text: str) -> set[str]:
    found = set()
    for raw in NUMBER.findall(ISO_DATE.sub(" ", text)):
        value = raw.replace(",", "")
        found.add(value)
        try:
            number = float(value)
            found.add(repr(number))
            if number.is_integer():
                found.add(str(int(number)))
        except ValueError:
            pass
    return found


def _significant(value: str) -> bool:
    """Small whole numbers ("two genomes", "step 3") are not checked."""
    try:
        number = float(value)
    except ValueError:
        return False
    return not (number.is_integer() and 0 <= number < 10)


def check(text: str, corpus: str, project: Project) -> list[str]:
    """Return a list of problems; an empty list means the answer passed."""
    problems = []
    unknown_ids = sorted(set(PROTEIN_ID.findall(text)) - set(PROTEIN_ID.findall(corpus)))
    if unknown_ids:
        problems.append("protein IDs not found in any tool result: " + ", ".join(unknown_ids))
    supported = _numbers(PROTEIN_ID.sub(" ", corpus))
    unsupported = sorted({raw.replace(",", "") for raw in NUMBER.findall(ISO_DATE.sub(" ", PROTEIN_ID.sub(" ", JOB_ID.sub(" ", text))))
                          if _significant(raw.replace(",", "")) and not (_numbers(raw) & supported)})
    if unsupported:
        problems.append("numbers not found in any tool result: " + ", ".join(unsupported[:10]))
    jobs = {j["id"]: j for j in project.jobs()}
    for job_id in sorted(set(JOB_ID.findall(text))):
        if job_id not in jobs:
            problems.append(f"job {job_id} does not exist")
    latest: dict[tuple[str, str], dict] = {}
    for job in jobs.values():
        for sample in [job.get("sample")] + list(job.get("samples") or []):
            if sample:
                latest[(sample, job["kind"])] = job
    samples = sorted({s for s, _ in latest}, key=len, reverse=True)
    for sentence in _sentences(text):
        if not COMPLETION.search(sentence) or HEDGE.search(sentence):
            continue
        lowered = sentence.lower()
        for job_id in JOB_ID.findall(sentence):
            state = jobs.get(job_id, {}).get("status")
            if state and state != "COMPLETE":
                problems.append(f"says job {job_id} finished, but its recorded status is {state}")
        kinds = [k for k, word in (("annotation", "annotat"), ("comparison", "compar"), ("ncbi_revision", "ncbi")) if word in lowered]
        named = [s for s in samples if re.search(rf"(?<![\w-]){re.escape(s)}(?![\w-])", sentence)]
        for sample in named:
            for kind in kinds or ["annotation"]:
                job = latest.get((sample, kind))
                if job and job.get("status") != "COMPLETE":
                    problems.append(f"says the {kind} of {sample} finished, but its latest {kind} job is {job.get('status')}")
                elif job is None and kinds:
                    problems.append(f"says the {kind} of {sample} finished, but there is no such job")
        if not named:
            for kind in kinds:
                if not any(j.get("kind") == kind and j.get("status") == "COMPLETE" for j in jobs.values()):
                    problems.append(f"claims a {kind} finished, but no {kind} job is COMPLETE")
    for match in NOVEL.finditer(text):
        if not NEGATION.search(text[max(0, match.start() - 60):match.start()]):
            problems.append("uses 'novel'; unknown function or no database match does not establish novelty")
            break
    return problems


def verified_fallback(results: list[dict]) -> str:
    """Shown instead of an answer that failed grounding twice: the tool results themselves."""
    if not results:
        return ("I could not verify my answer against PhageMine's outputs, so I have withheld it. "
                "Ask again and I will look the results up.")
    lines = ["I could not verify my written answer against PhageMine's outputs, so here is exactly what the tools returned:"]
    for item in results:
        result = item["result"]
        if not result.get("ok"):
            lines.append(f"\n**{item['name']}** reported an error: {result.get('error')}")
            continue
        lines.append(f"\n**{item['name']}**")
        for key, value in result.items():
            if key == "ok":
                continue
            if isinstance(value, list) and value and isinstance(value[0], dict):
                lines.append(f"- {key}: {len(value)} rows")
                columns = list(value[0])[:6]
                lines.append("\n| " + " | ".join(columns) + " |\n| " + " | ".join("---" for _ in columns) + " |")
                for row in value[:15]:
                    lines.append("| " + " | ".join(str(row.get(c, "")) for c in columns) + " |")
            else:
                rendered = json.dumps(value, default=str) if isinstance(value, (dict, list)) else str(value)
                lines.append(f"- {key}: {rendered[:300]}")
    return "\n".join(lines)
