"""Deterministic checks that an answer only states what tools returned.

The check is deliberately narrow and strict: protein IDs, job IDs, completion
claims, and the unknown-means-novel error. It does not judge prose quality.
"""
from __future__ import annotations

import re

from .store import Project

PROTEIN_ID = re.compile(r"\bPM_\d{6}\b")
JOB_ID = re.compile(r"\b(?:annotation|comparison|ncbirevision)-[0-9a-f]{8}\b")
COMPLETION = re.compile(r"\b(?:completed?|finished|succeeded|is done|are done|has been (?:annotated|compared))\b", re.I)
NOVEL = re.compile(r"\bnovel\b", re.I)
NEGATION = re.compile(r"(?:not|n't|no|never|does not mean|doesn't mean|without)\W+(?:\w+\W+){0,4}$", re.I)


def _sentences(text: str) -> list[str]:
    return [s for s in re.split(r"(?<=[.!?\n])\s+", text) if s.strip()]


def check(text: str, corpus: str, project: Project) -> list[str]:
    """Return a list of problems; an empty list means the answer passed."""
    problems = []
    unknown_ids = sorted(set(PROTEIN_ID.findall(text)) - set(PROTEIN_ID.findall(corpus)))
    if unknown_ids:
        problems.append("protein IDs not found in any tool result: " + ", ".join(unknown_ids))
    jobs = {j["id"]: j for j in project.jobs()}
    for job_id in sorted(set(JOB_ID.findall(text))):
        if job_id not in jobs:
            problems.append(f"job {job_id} does not exist")
    for sentence in _sentences(text):
        if not COMPLETION.search(sentence):
            continue
        for job_id in JOB_ID.findall(sentence):
            state = jobs.get(job_id, {}).get("status")
            if state and state != "COMPLETE" and not re.search(r"\bnot\b|\bfailed\b|\brunning\b|\bwhen\b|\bonce\b|\bafter\b", sentence, re.I):
                problems.append(f"sentence says job {job_id} finished, but its recorded status is {state}")
        lowered = sentence.lower()
        for kind, word in (("annotation", "annotat"), ("comparison", "compar")):
            if word in lowered and not any(j.get("kind") == kind and j.get("status") == "COMPLETE" for j in jobs.values()):
                if not re.search(r"\b(?:will|once|when|after|until|not yet|still|running|queued|waiting|failed)\b", lowered):
                    problems.append(f"claims a {kind} finished, but no {kind} job is COMPLETE")
    for match in NOVEL.finditer(text):
        before = text[max(0, match.start() - 60):match.start()]
        if not NEGATION.search(before):
            problems.append("uses 'novel'; unknown function or no database match does not establish novelty")
            break
    return problems
