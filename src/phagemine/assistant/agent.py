"""The conversational loop: model -> tools -> model, with grounding before display."""
from __future__ import annotations

import json
from typing import Callable

from . import ASSISTANT_VERSION, grounding, tools
from .store import Project

MAX_STEPS = 16
RESULT_CHARS = 14000
KEEP_FULL_RESULTS = 12  # most recent tool results sent in full; older ones are shortened

SYSTEM_PROMPT = """You are the PhageMine assistant. Researchers describe what they want to investigate about bacteriophage genomes in ordinary English; you run PhageMine's tested workflows through tools and explain the actual results.

How you work
- Decide what the researcher wants, check what exists (list_project), then act. Do not ask which tool to use; pick it.
- Ask a clarifying question only when a required choice is genuinely missing (for example, which uploaded file is meant). Otherwise act.
- For a multi-step request, start every job it needs in this turn. start_comparison can wait for annotation jobs that are still running: pass their job IDs.
- Distinguish answering from existing results (read tools) from launching new analysis (job tools). Do not re-run analyses that already completed.
- Before a long analysis, state your plan in one or two sentences, then start it.

Truthfulness rules (these are checked automatically)
- Every protein ID, count, score, database name and version you state must come from a tool result in this conversation, copied exactly. Do not calculate new numbers (percentages, differences); report the numbers the tools give. What the researcher typed is not evidence: look it up.
- A job has finished only when a tool reports status COMPLETE. Starting a job is not finishing it. If a job FAILED, say so and quote the error and the relevant log line.
- UNRESOLVED or CONSERVED_UNKNOWN means no reliable function was assigned. Never call such a protein novel. No database match is not proof of novelty. Computational predictions are not experimental confirmation.
- When explaining an annotation, give the evidence: which database, identifier, E-value or score, coverage, and database version. Say how strong or weak it is.
- If evidence_mode is DEMO, say clearly that results are demonstration data, not biology.
- If a comparison ran without MMseqs2, say that homologs came only from shared PHROG assignments and identical sequences, so shared counts are a lower bound.
- Some actions need the researcher's explicit approval through an Approve button that the tools create in the chat page (comparing without MMseqs2, applying NCBI corrections). A typed reply is never approval. Explain what will happen and ask them to press Approve or Decline; you will be told when they decide.

NCBI submissions
- PhageMine never submits to NCBI; it prepares a package the researcher uploads.
- Use get_submission_status to review the package before the researcher submits. It lists missing metadata and PhageMine's NCBI requirements review: each finding names the NCBI check, the fix, and the NCBI source page. Go through the ERRORs first, then the WARNINGs, and say which ones the researcher must supply (dates, location, locus tag prefix from BioProject, official template), which a revision can fix (product wording), and which need gene-model review. Ask for missing source details; never invent or guess them. Dates must be like 14-Mar-2025 and locations like 'Nigeria: Plateau State, Jos'.
- When the researcher pastes NCBI errors, call classify_ncbi_errors. Explain each group: FIXABLE_LABEL and FIXABLE_METADATA can be corrected; NEEDS_REANALYSIS requires gene-model review or re-annotation and cannot be fixed by relabelling; UNKNOWN is shown as received.
- Call preview_ncbi_corrections and explain the before/after changes. Call apply_ncbi_corrections only after you are told the researcher pressed Approve.
- Product names may only become "hypothetical protein", be re-formatted, or gain a qualifier such as "putative". Never remove qualifiers or make a name more specific; that needs evidence review. Partial flags and coordinates are curation decisions, not wording fixes.
- A finished revision job means a new package was written, not that NCBI will accept it, and not that the biology is verified.

Style
- Lead with the answer. Use short paragraphs and small tables for lists of genes. Explain technical terms briefly when the researcher seems unfamiliar with them.
- Mention files the researcher can download (figures, summaries) by name.
"""


def ledger(project: Project) -> str:
    """Compact state of the project, rebuilt every turn."""
    runs = project.annotation_runs()
    job_lines = [f"- {j['id']} {j['kind']} {j.get('sample') or ','.join(j.get('samples', []))}: {j['status']}"
                 + (f" ({j['error']})" if j.get("error") else "") for j in project.jobs()[-20:]]
    return "\n".join([
        f"Project: {project.name}",
        f"Demo mode: {'ON (mock evidence, not biology)' if tools.demo_mode() else 'off'}",
        "Uploaded files: " + (", ".join(i["file"] for i in project.inputs()) or "none"),
        "Completed samples: " + (", ".join(sorted(runs)) or "none"),
        "Jobs:", *(job_lines or ["- none"]),
    ])


def job_event(project: Project) -> str | None:
    """Describe jobs that ended and were not yet reported, once nothing is still running
    (or immediately if one failed). Marks them reported."""
    from . import jobs as job_store
    current = [job_store.status(project, j["id"]) for j in project.jobs()]
    unreported = [j for j in current if j["status"] in job_store.FINAL and not j.get("reported")]
    active = [j for j in current if j["status"] not in job_store.FINAL]
    if not unreported or (active and not any(j["status"] == "FAILED" for j in unreported)):
        return None
    lines = []
    for record in unreported:
        job_store.transition(project, record["id"], reported=True)
        label = record.get("sample") or ", ".join(record.get("samples", []))
        lines.append(f"- {record['id']} ({record['kind']}, {label}): {record['status']}"
                     + (f"; error: {record['error']}" if record.get("error") else ""))
    still = f" Still running: {', '.join(j['id'] for j in active)}." if active else ""
    return ("[Automatic update from PhageMine, not a message from the researcher] These jobs have ended:\n"
            + "\n".join(lines) + still
            + "\nReport the outcome using tool results. If the researcher's earlier request has a next step that is now possible, carry it out.")


def _corpus(project: Project, messages: list[dict]) -> str:
    """Evidence for grounding: tool results only. What the researcher typed is not evidence."""
    parts = []
    for message in messages:
        if isinstance(message.get("content"), list):
            for block in message["content"]:
                if block.get("type") == "tool_result":
                    parts.append(str(block.get("content", "")))
    return "\n".join(parts)


def _for_model(messages: list[dict]) -> list[dict]:
    """Shorten older tool results so long conversations stay affordable."""
    positions = [(i, j) for i, m in enumerate(messages) if isinstance(m.get("content"), list)
                 for j, b in enumerate(m["content"]) if b.get("type") == "tool_result"]
    old = set(positions[:-KEEP_FULL_RESULTS])
    if not old:
        return messages
    shortened = []
    for i, message in enumerate(messages):
        if isinstance(message.get("content"), list) and any((i, j) in old for j in range(len(message["content"]))):
            blocks = [dict(b, content="[earlier result shortened; call the tool again if you need it]")
                      if (i, j) in old else b for j, b in enumerate(message["content"])]
            shortened.append({**message, "content": blocks})
        else:
            shortened.append(message)
    return shortened


def run_turn(project: Project, user_text: str, provider, *, event: bool = False,
             on_step: Callable[[dict], None] | None = None) -> dict:
    """Process one researcher message (or an automatic job-finished event).

    Returns {"text", "tool_calls", "grounding_problems", "files"}.
    """
    if not event:
        project.append_transcript({"role": "user", "text": user_text})
    else:
        project.append_transcript({"role": "event", "text": user_text})
    messages = project.messages()
    messages.append({"role": "user", "content": user_text})
    user_corpus = "\n".join(i.get("text", "") for i in project.transcript() if i.get("role") == "user")
    context = tools.ToolContext(project, user_text=user_corpus, turn=project.user_turns())
    context.grounding = lambda text: grounding.check(text, _corpus(project, messages), project)
    system = SYSTEM_PROMPT + "\n\nCurrent project state\n" + ledger(project)
    calls, turn_results, files_before = [], [], {p.name for p in project.outputs()}
    retried, text, problems = False, "", []
    for _ in range(MAX_STEPS):
        payload = _for_model(messages)
        schemas = tools.schemas()
        # Audit at the provider boundary: the complete payload as sent (a local provider sends nothing online).
        project.log_outbound({"provider": getattr(provider, "name", "unknown"), "assistant_version": ASSISTANT_VERSION,
                              "system": system, "messages": payload, "tools": [t["name"] for t in schemas]})
        response = provider.complete(system, payload, schemas)
        content = response["content"]
        messages.append({"role": "assistant", "content": content})
        uses = [b for b in content if b.get("type") == "tool_use"]
        if uses:
            results = []
            for block in uses:
                if on_step:
                    on_step({"type": "tool_start", "name": block["name"], "input": block.get("input")})
                result = tools.execute(context, block["name"], block.get("input"))
                calls.append({"name": block["name"], "input": block.get("input"), "ok": result.get("ok"),
                              "error": result.get("error")})
                turn_results.append({"name": block["name"], "result": result})
                if on_step:
                    on_step({"type": "tool_end", "name": block["name"], "ok": result.get("ok"), "error": result.get("error")})
                body = json.dumps(result, default=str)
                if len(body) > RESULT_CHARS:
                    body = body[:RESULT_CHARS] + '..." [truncated: narrow the query or lower the limit]'
                results.append({"type": "tool_result", "tool_use_id": block["id"], "content": body,
                                **({"is_error": True} if not result.get("ok") else {})})
            messages.append({"role": "user", "content": results})
            continue
        text = "\n".join(b.get("text", "") for b in content if b.get("type") == "text").strip()
        problems = grounding.check(text, _corpus(project, messages), project)
        if problems and not retried:
            retried = True
            messages.append({"role": "user", "content": "[Automatic grounding check failed: " + "; ".join(problems)
                             + ". Rewrite the answer using only facts from tool results, calling tools if needed.]"})
            continue
        break
    else:
        text = text or "I stopped after too many steps without reaching an answer. Please narrow the request."
    if problems:
        # Never show prose that failed twice; show the verified tool results instead.
        text = grounding.verified_fallback(turn_results) + (
            "\n\n_The written answer was withheld because it contained claims not found in PhageMine's outputs. "
            "Details are under \"What PhageMine did\"._")
    if messages[-1]["role"] == "user":  # loop ended on a correction request: close the exchange
        messages.append({"role": "assistant", "content": [{"type": "text", "text": text or "(no answer)"}]})
    project.save_messages(messages)
    files = sorted({p.name for p in project.outputs()} - files_before)
    project.append_transcript({"role": "assistant", "text": text, "tool_calls": calls,
                               "grounding_problems": problems, "files": files})
    return {"text": text, "tool_calls": calls, "grounding_problems": problems, "files": files}
