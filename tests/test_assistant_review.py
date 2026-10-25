"""Counterexamples from external review of PR #39. Each asserts the correct behaviour."""
from __future__ import annotations

import json
import os
import stat
import time
from pathlib import Path

import pytest

from phagemine.assistant import agent, approvals, grounding, jobs, ncbi, runner, tools
from phagemine.assistant.llm import ScriptedProvider
from phagemine.assistant.store import Project

DEMO = Path(__file__).resolve().parents[1] / "examples" / "demo_phage.fasta"


@pytest.fixture(autouse=True)
def demo_env(tmp_path, monkeypatch):
    monkeypatch.setenv("PHAGEMINE_ASSISTANT_DEMO", "1")
    monkeypatch.setenv("PHAGEMINE_PROJECTS", str(tmp_path / "projects"))
    monkeypatch.setattr(tools.shutil, "which", lambda name: None)


def rotated() -> bytes:
    seq = "".join(DEMO.read_text().splitlines()[1:])
    return f">demo_rotated\n{seq[len(seq)//3:]}{seq[:len(seq)//3]}\n".encode()


def ctx(p, turn=1, text=""):
    return tools.ToolContext(p, user_text=text, turn=turn)


def call(c, name, **arguments):
    return tools.execute(c, name, arguments)


def wait(p, job_id, timeout=120):
    deadline = time.time() + timeout
    while time.time() < deadline:
        record = jobs.status(p, job_id)
        if record["status"] in jobs.FINAL:
            return record
        time.sleep(0.5)
    raise AssertionError(jobs.status(p, job_id))


def consent_exact_only(p):
    request = approvals.request(p, "exact_only_comparison", {"reason": "MMseqs2 unavailable"}, "Compare without MMseqs2")
    approvals.approve(p, request["id"])


@pytest.fixture
def compared():
    p = Project.create("review")
    p.add_input("demo_phage.fasta", DEMO.read_bytes())
    p.add_input("demo_rotated.fasta", rotated())
    started = call(ctx(p), "start_annotation", files=["demo_phage.fasta", "demo_rotated.fasta"])
    ids = [s["job_id"] for s in started["started"]]
    consent_exact_only(p)
    comparison = call(ctx(p), "start_comparison", annotation_job_ids=ids, allow_exact_only=True)
    assert comparison["ok"], comparison
    assert wait(p, comparison["job_id"])["status"] == "COMPLETE"
    return p, comparison["job_id"]


# 1. Completion requires validated outputs --------------------------------

def test_malformed_annotation_output_fails_validation(compared):
    p, _ = compared
    record = next(j for j in p.jobs() if j["kind"] == "annotation")
    table = p.root / record["output"] / "annotation.tsv"
    lines = table.read_text().splitlines()
    table.chmod(stat.S_IWUSR | stat.S_IRUSR)
    table.write_text("\n".join(lines[:-1]) + "\n")
    with pytest.raises(ValueError):
        runner.validate_outputs(p, record)


def test_invalid_comparison_json_fails_validation(compared):
    p, comparison_id = compared
    record = p.job(comparison_id)
    (p.root / record["directions"][0]["output"] / "comparative_context.json").write_text("not json")
    with pytest.raises(ValueError):
        runner.validate_outputs(p, record)


def test_worker_dead_while_queued_is_failed():
    p = Project.create("queued")
    record = jobs.create_job(p, "annotation", [["phagemine", "--version"]], output="runs/x", launch=False, meta={"sample": "x"})
    stored = p.job(record["id"]); stored.update(worker_pid=999999, created="2020-01-01T00:00:00+00:00"); p.save_job(stored)
    assert jobs.status(p, record["id"])["status"] == "FAILED"


def test_final_states_are_never_overwritten():
    p = Project.create("race")
    record = jobs.create_job(p, "annotation", [["phagemine", "--version"]], output="runs/x", launch=False, meta={"sample": "x"})
    jobs.cancel(p, record["id"])
    assert jobs.transition(p, record["id"], status="COMPLETE")["status"] == "CANCELLED"


# 2. Grounding -------------------------------------------------------------

def test_completion_claim_is_bound_to_the_named_sample(compared):
    p, _ = compared
    failed = jobs.create_job(p, "annotation", [["phagemine", "annotate", "inputs/none.fasta", "--output", "runs/{job_id}"]],
                             output="runs/{job_id}", meta={"sample": "phageB"})
    wait(p, failed["id"])
    problems = grounding.check("Annotation of phageB completed successfully.", "", p)
    assert any("phageB" in problem for problem in problems)


def test_unsupported_numbers_are_flagged(compared):
    p, _ = compared
    corpus = json.dumps({"proteins_with_homolog_in_other": {"demo_phage": 4}})
    assert grounding.check("There are 4 shared proteins.", corpus, p) == []
    assert any("999999" in problem for problem in grounding.check("There are 999999 shared genes.", corpus, p))


def test_user_text_is_not_evidence():
    p = Project.create("user text")
    p.append_transcript({"role": "user", "text": "Is PM_000123 a capsid?"})
    provider = ScriptedProvider([[{"type": "text", "text": "Yes, PM_000123 is a capsid."}]] * 2)
    result = agent.run_turn(p, "Is PM_000123 a capsid?", provider)
    assert result["grounding_problems"]


def test_rejected_answer_is_not_displayed_after_retry():
    p = Project.create("fallback")
    provider = ScriptedProvider([[{"type": "text", "text": "There are 999999 shared genes."}]] * 2)
    result = agent.run_turn(p, "How many genes are shared?", provider)
    assert "999999" not in result["text"]
    assert "could not verify" in result["text"].lower()


# 3. NCBI ------------------------------------------------------------------

def test_later_message_is_not_approval(compared):
    p, _ = compared
    preview = call(ctx(p, 2, "host is Escherichia coli"), "preview_ncbi_corrections", sample="demo_phage",
                   metadata={"host": "Escherichia coli"})
    assert preview["ok"], preview
    refused = call(ctx(p, 3, "host is Escherichia coli\nDo not apply the corrections."), "apply_ncbi_corrections",
                   preview_id=preview["preview_id"])
    assert refused["ok"] is False and "approv" in refused["error"].lower()
    approvals.approve(p, preview["approval_id"])
    applied = call(ctx(p, 3, "host is Escherichia coli"), "apply_ncbi_corrections", preview_id=preview["preview_id"])
    assert applied["ok"], applied


def test_approval_is_bound_to_the_exact_corrections(compared):
    p, _ = compared
    preview = call(ctx(p, 2, "host is Escherichia coli"), "preview_ncbi_corrections", sample="demo_phage",
                   metadata={"host": "Escherichia coli"})
    approvals.approve(p, preview["approval_id"])
    state = p.read_state(f"ncbi/{preview['preview_id']}.json", {})
    state["corrections"]["metadata"]["host"] = "Salmonella"
    p.write_state(f"ncbi/{preview['preview_id']}.json", state)
    refused = call(ctx(p, 3), "apply_ncbi_corrections", preview_id=preview["preview_id"])
    assert refused["ok"] is False and "changed" in refused["error"]


def test_uncertainty_qualifiers_are_preserved():
    assert ncbi.check_product("DNA polymerase", "putative DNA polymerase", []) is not None
    assert ncbi.check_product("hypothetical protein", "putative DNA polymerase", []) is None
    assert ncbi.check_product("putative DNA polymerase", "Putative DNA polymerase.", []) is None


def test_free_text_evidence_cannot_authorise_a_name():
    statement = "Evidence is insufficient to assign DNA polymerase function."
    assert ncbi.check_product("DNA polymerase", "hypothetical protein", [statement]) is not None


def test_partial_flags_are_not_routine_corrections(compared):
    p, _ = compared
    result = call(ctx(p, 2), "preview_ncbi_corrections", sample="demo_phage", annotations={"PM_000001": {"partial": True}})
    assert result["ok"] is False and "curation" in result["error"]


def test_unrecognised_validator_lines_are_kept():
    result = ncbi.classify("Line one about something\nAnother remark from the curator")
    assert result["total_messages"] == 2
    assert {m["outcome"] for m in result["messages"]} == {ncbi.UNKNOWN}


def test_successive_corrections_merge_field_by_field(compared):
    p, _ = compared
    first = call(ctx(p, 2), "preview_ncbi_corrections", sample="demo_phage",
                 annotations={"PM_000001": {"product": "hypothetical protein", "reason": "a"}})
    approvals.approve(p, first["approval_id"])
    wait(p, call(ctx(p, 3), "apply_ncbi_corrections", preview_id=first["preview_id"])["job_id"])
    second = call(ctx(p, 4), "preview_ncbi_corrections", sample="demo_phage",
                  annotations={"PM_000001": {"note": "checked", "reason": "b"}})
    merged = p.read_state(f"ncbi/{second['preview_id']}.json", {})["corrections"]["annotations"]["PM_000001"]
    assert merged["product"] == "hypothetical protein" and merged["note"] == "checked"


def test_existing_submission_metadata_is_preserved(compared):
    p, _ = compared
    package = p.run_dir("demo_phage") / "genbank_submission" / "submission_metadata.json"
    data = json.loads(package.read_text()); data["organism"] = "Escherichia phage T7-like"
    package.write_text(json.dumps(data))
    preview = call(ctx(p, 2, "host is Escherichia coli"), "preview_ncbi_corrections", sample="demo_phage",
                   metadata={"host": "Escherichia coli"})
    merged = p.read_state(f"ncbi/{preview['preview_id']}.json", {})["corrections"]["metadata"]
    assert merged["organism"] == "Escherichia phage T7-like" and merged["host"] == "Escherichia coli"


# 4. Comparison ------------------------------------------------------------

def test_mmseqs_reported_from_execution_not_presence(tmp_path, monkeypatch):
    fake = tmp_path / "mmseqs"
    fake.write_text("#!/bin/sh\necho broken >&2\nexit 1\n"); fake.chmod(0o755)
    monkeypatch.setattr(tools.shutil, "which", lambda name: str(fake))
    p = Project.create("broken mmseqs")
    p.add_input("demo_phage.fasta", DEMO.read_bytes()); p.add_input("demo_rotated.fasta", rotated())
    ids = [s["job_id"] for s in call(ctx(p), "start_annotation", files=["demo_phage.fasta", "demo_rotated.fasta"])["started"]]
    comparison = call(ctx(p), "start_comparison", annotation_job_ids=ids)
    wait(p, comparison["job_id"])
    shared = call(ctx(p), "get_shared_genes")
    assert shared["mmseqs2_used"] is False and "warning" in shared


def test_pair_table_agrees_with_counts(compared):
    p, comparison_id = compared
    record = p.job(comparison_id)
    forward = p.root / record["directions"][0]["output"] / "ortholog_links.tsv"
    header = forward.read_text().splitlines()[0]
    forward.write_text(header + "\n")  # links now exist only in the reverse direction
    shared = call(ctx(p), "get_shared_genes")
    a, b = shared["sample_a"], shared["sample_b"]
    assert shared["proteins_with_homolog_in_other"][a] > 0
    assert shared["pairs_total"] > 0
    assert all(pair["support"] == f"{b}->{a}" for pair in shared["pairs"])


def test_comparison_uses_the_runs_it_compared(compared):
    p, comparison_id = compared
    rerun = call(ctx(p), "start_annotation", files=["demo_phage.fasta"], rerun=True)
    wait(p, rerun["started"][0]["job_id"])
    shared = call(ctx(p), "get_shared_genes")
    assert shared["source_runs"]["demo_phage"] == p.job(comparison_id)["run_dirs"]["demo_phage"]


def test_exact_only_requires_recorded_consent():
    p = Project.create("consent")
    p.add_input("demo_phage.fasta", DEMO.read_bytes()); p.add_input("demo_rotated.fasta", rotated())
    ids = [s["job_id"] for s in call(ctx(p), "start_annotation", files=["demo_phage.fasta", "demo_rotated.fasta"])["started"]]
    refused = call(ctx(p), "start_comparison", annotation_job_ids=ids, allow_exact_only=True)
    assert refused["ok"] is False and "approval" in refused["error"].lower()
    assert approvals.pending(p)


# 5. Cross-cutting ---------------------------------------------------------

def test_argument_types_and_bounds_are_enforced():
    p = Project.create("schema")
    p.add_input("demo_phage.fasta", DEMO.read_bytes())
    assert call(ctx(p), "start_annotation", files=["demo_phage.fasta"], threads=10000)["ok"] is False
    assert call(ctx(p), "start_annotation", files="demo_phage.fasta")["ok"] is False
    assert call(ctx(p), "query_proteins", functional_states=["MADE_UP"])["ok"] is False


def test_audit_log_records_the_full_payload():
    p = Project.create("audit")
    agent.run_turn(p, "hello", ScriptedProvider([[{"type": "text", "text": "Hello."}]]))
    entry = json.loads((p.root / "outbound.jsonl").read_text().splitlines()[-1])
    assert entry["system"].startswith("You are the PhageMine assistant")
    assert entry["messages"][-1]["content"] == "hello"
    assert "list_project" in entry["tools"]
