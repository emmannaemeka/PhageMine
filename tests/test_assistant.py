"""Assistant tests: tools, jobs, grounding, NCBI guards and the agent loop.

Analyses run in demo mode (fixture gene caller, mock evidence), so these tests
check orchestration and honesty, not biology.
"""
from __future__ import annotations

import json
import os
import time
from pathlib import Path

import pytest

from phagemine.assistant import agent, grounding, jobs, ncbi, tools
from phagemine.assistant.llm import ScriptedProvider
from phagemine.assistant.store import Project

DEMO = Path(__file__).resolve().parents[1] / "examples" / "demo_phage.fasta"


def rotated_demo() -> bytes:
    lines = DEMO.read_text().splitlines()
    seq = "".join(lines[1:])
    third = len(seq) // 3
    return f">demo_rotated\n{seq[third:]}{seq[:third]}\n".encode()


@pytest.fixture(autouse=True)
def demo_env(tmp_path, monkeypatch):
    monkeypatch.setenv("PHAGEMINE_ASSISTANT_DEMO", "1")
    monkeypatch.setenv("PHAGEMINE_PROJECTS", str(tmp_path / "projects"))
    monkeypatch.setattr(tools.shutil, "which", lambda name: None)  # MMseqs2 absent, deterministically


@pytest.fixture
def project():
    p = Project.create("test project")
    p.add_input("demo_phage.fasta", DEMO.read_bytes())
    p.add_input("demo_rotated.fasta", rotated_demo())
    return p


def ctx(p, turn=1, text=""):
    return tools.ToolContext(p, user_text=text, turn=turn)


def call(context, name, **arguments):
    return tools.execute(context, name, arguments)


def wait(p, job_id, timeout=120):
    deadline = time.time() + timeout
    while time.time() < deadline:
        record = jobs.status(p, job_id)
        if record["status"] in jobs.FINAL:
            return record
        time.sleep(0.5)
    raise AssertionError(f"job {job_id} did not finish: {jobs.status(p, job_id)}")


@pytest.fixture
def compared(project):
    c = ctx(project)
    started = call(c, "start_annotation", files=["demo_phage.fasta", "demo_rotated.fasta"])
    ids = [s["job_id"] for s in started["started"]]
    comparison = call(c, "start_comparison", annotation_job_ids=ids, allow_exact_only=True)
    assert comparison["ok"], comparison
    assert wait(project, comparison["job_id"])["status"] == "COMPLETE"
    return project, ids, comparison["job_id"]


# ------------------------------------------------------------------ store ----

def test_inputs_are_copied_read_only_and_paths_are_contained(project):
    path = project.input_path("demo_phage.fasta")
    assert not os.access(path, os.W_OK) or os.geteuid() == 0
    with pytest.raises(ValueError, match="different content"):
        project.add_input("demo_phage.fasta", b">x\nACGT\n")
    with pytest.raises(ValueError, match="outside"):
        project.contained("../../etc/passwd")
    with pytest.raises(ValueError):
        project.input_path("../project.json")


def test_projects_are_isolated(compared):
    other = Project.create("other")
    with pytest.raises(ValueError, match="No completed annotation"):
        other.run_dir("demo_phage")
    result = call(ctx(other), "get_run_summary", sample="demo_phage")
    assert result["ok"] is False


# ------------------------------------------------------------------- jobs ----

def test_annotation_job_completes_only_with_validated_outputs(project):
    started = call(ctx(project), "start_annotation", files=["demo_phage.fasta"])
    record = wait(project, started["started"][0]["job_id"])
    assert record["status"] == "COMPLETE" and record["exit_codes"] == [0]
    assert record["pipeline_version"]


def test_nonzero_exit_is_failed(project):
    record = jobs.create_job(project, "annotation", [["phagemine", "annotate", "inputs/missing.fasta", "--output", "runs/{job_id}"]],
                             output="runs/{job_id}", meta={"sample": "missing"})
    final = wait(project, record["id"])
    assert final["status"] == "FAILED" and "exited with code" in final["error"]


def test_exit_zero_without_artifacts_is_failed(project):
    record = jobs.create_job(project, "annotation", [["phagemine", "--version"]], output="runs/{job_id}", meta={"sample": "x"})
    final = wait(project, record["id"])
    assert final["status"] == "FAILED" and "did not validate" in final["error"]


def test_failed_prerequisite_fails_dependent_job(project):
    bad = jobs.create_job(project, "annotation", [["phagemine", "annotate", "inputs/none.fasta", "--output", "runs/{job_id}"]],
                          output="runs/{job_id}", meta={"sample": "none"})
    dependent = jobs.create_job(project, "comparison", [["phagemine", "--version"]], output="runs/{job_id}",
                                depends_on=[bad["id"]], meta={"samples": [], "directions": []})
    final = wait(project, dependent["id"])
    assert final["status"] == "FAILED" and "prerequisite" in final["error"]


def test_dead_worker_is_reported_as_failed(project):
    record = jobs.create_job(project, "annotation", [["phagemine", "--version"]], output="runs/x", launch=False, meta={"sample": "x"})
    stored = project.job(record["id"]); stored.update(status="RUNNING", worker_pid=999999); project.save_job(stored)
    assert jobs.status(project, record["id"])["status"] == "FAILED"


def test_inputs_unchanged_after_jobs(compared):
    p, _, _ = compared
    before = {i["file"]: i["sha256"] for i in p.inputs()}
    assert before == {i["file"]: i["sha256"] for i in Project(p.root).inputs()}
    for job in p.jobs():
        assert all("inputs/" not in part or part.startswith("inputs/") for step in job["steps"] for part in step)


# ------------------------------------------------------------------ tools ----

def test_tool_arguments_are_validated(project):
    assert "unexpected" in call(ctx(project), "inspect_fasta", file="demo_phage.fasta", path="/etc")["error"]
    assert "missing" in tools.execute(ctx(project), "inspect_fasta", {})["error"]
    assert call(ctx(project), "no_such_tool")["ok"] is False


def test_comparison_refuses_exact_only_without_consent(project):
    started = call(ctx(project), "start_annotation", files=["demo_phage.fasta", "demo_rotated.fasta"])
    ids = [s["job_id"] for s in started["started"]]
    refused = call(ctx(project), "start_comparison", annotation_job_ids=ids)
    assert refused["ok"] is False and "MMseqs2" in refused["error"]


def test_shared_genes_are_counted_from_both_sides_with_method(compared):
    p, _, comparison_id = compared
    shared = call(ctx(p), "get_shared_genes")
    assert shared["comparison_job_id"] == comparison_id
    assert set(shared["proteins_with_homolog_in_other"]) == {"demo_phage", "demo_rotated"}
    assert all(v > 0 for v in shared["proteins_with_homolog_in_other"].values())
    assert shared["mmseqs2_used"] is False and "lower bound" in shared["warning"]
    assert set(shared["detection_methods"]) == {"EXACT_SEQUENCE_MATCH"}


def test_unknown_function_query_and_evidence(compared):
    p, _, _ = compared
    unknown = call(ctx(p), "query_proteins", unknown_function_only=True)
    assert unknown["total_matches"] > 0
    assert all(r["functional_state"] in tools.UNKNOWN_STATES for r in unknown["proteins"])
    first = unknown["proteins"][0]
    evidence = call(ctx(p), "get_protein_evidence", sample=first["sample"], protein_id=first["protein_id"])
    assert evidence["evidence_mode"].startswith("DEMO")
    assert "protein_sequence" not in evidence


def test_comparison_figure_has_source_table(compared):
    p, _, _ = compared
    result = call(ctx(p), "make_comparison_figure")
    assert result["ok"], result
    table = p.root / "figures" / next(f for f in result["files"] if f.endswith(".tsv"))
    rows = table.read_text().splitlines()[1:]
    assert sum(r.endswith("True") for r in rows) == sum(result["shared_counts"].values())


def test_summary_refused_when_not_grounded(compared):
    p, _, _ = compared
    c = ctx(p)
    c.grounding = lambda text: grounding.check(text, "", p)
    refused = call(c, "write_summary", title="Summary", markdown="PM_999999 is a novel capsid.")
    assert refused["ok"] is False and "not grounded" in refused["error"]
    c.grounding = lambda text: []
    saved = call(c, "write_summary", title="Summary", markdown="Two genomes were compared.")
    assert saved["ok"] and (p.root / "summaries" / saved["file"]).read_text().count("DEMO MODE") == 1


# ------------------------------------------------------------------- NCBI ----

def test_ncbi_errors_are_routed_to_outcomes():
    result = ncbi.classify("ERROR: SEQ_FEAT.InternalStop in PM_000004\nWARNING: missing collection-date\n"
                           "Suspect product name in PM_000002\nERROR: something unheard of")
    outcomes = [m["outcome"] for m in result["messages"]]
    assert outcomes == [ncbi.NEEDS_REANALYSIS, ncbi.FIXABLE_METADATA, ncbi.FIXABLE_LABEL, ncbi.UNKNOWN]
    assert result["messages"][0]["protein_ids"] == ["PM_000004"]


def test_product_names_never_become_more_specific():
    assert ncbi.check_product("hypothetical protein", "putative tail fiber", []) is None
    assert ncbi.check_product("tail fiber", "putative tail fiber protein", []) is None
    assert ncbi.check_product("major capsid protein", "hypothetical protein", ["unrelated statement"])
    assert ncbi.check_product("terminase large subunit", "hypothetical protein", ["PHROG 1: terminase large subunit"]) is None


def test_ncbi_corrections_need_typed_metadata_and_later_approval(compared):
    p, _, _ = compared
    typed = "The host is Escherichia coli K-12 and the collection date is 2025-03-14."
    invented = call(ctx(p, 2, typed), "preview_ncbi_corrections", sample="demo_phage", metadata={"collection_date": "2024-01-01"})
    assert invented["ok"] is False and "not provided by the researcher" in invented["error"]
    coords = call(ctx(p, 2, typed), "preview_ncbi_corrections", sample="demo_phage", annotations={"PM_000001": {"start": 5}})
    assert coords["ok"] is False and "cannot be revised" in coords["error"]
    preview = call(ctx(p, 2, typed), "preview_ncbi_corrections", sample="demo_phage",
                   metadata={"host": "Escherichia coli K-12", "collection_date": "2025-03-14"},
                   annotations={"PM_000001": {"product": "hypothetical protein", "reason": "NCBI"}})
    assert preview["ok"], preview
    same_turn = call(ctx(p, 2, typed), "apply_ncbi_corrections", preview_id=preview["preview_id"])
    assert same_turn["ok"] is False and "not yet seen" in same_turn["error"]
    applied = call(ctx(p, 3, typed), "apply_ncbi_corrections", preview_id=preview["preview_id"])
    assert wait(p, applied["job_id"])["status"] == "COMPLETE"
    status = call(ctx(p, 3, typed), "get_submission_status", sample="demo_phage")
    assert status["package_version"].startswith("revision")
    assert "host" not in status["missing_required_metadata"]
    assert status["corrections_applied"]["metadata"]["host"] == "Escherichia coli K-12"
    again = call(ctx(p, 4, typed), "apply_ncbi_corrections", preview_id=preview["preview_id"])
    assert again["ok"] is False and "already applied" in again["error"]


# -------------------------------------------------------------- grounding ----

def test_grounding_flags_invented_ids_false_completion_and_novelty(project):
    corpus = json.dumps({"proteins": [{"protein_id": "PM_000001"}]})
    assert grounding.check("PM_000001 has no assigned function.", corpus, project) == []
    assert any("PM_000099" in p for p in grounding.check("PM_000099 is shared.", corpus, project))
    assert any("annotation finished" in p for p in grounding.check("The annotation is completed.", corpus, project))
    assert grounding.check("I will tell you once the annotation is completed.", corpus, project) == []
    assert any("novel" in p for p in grounding.check("PM_000001 is a novel protein.", corpus, project))
    assert grounding.check("Unknown function does not mean novel.", corpus, project) == []


# ------------------------------------------------------------------ agent ----

def tool_use(name, tool_id, **arguments):
    return {"type": "tool_use", "id": tool_id, "name": name, "input": arguments}


def text(value):
    return {"type": "text", "text": value}


def test_agent_runs_tools_and_reports_only_starting(project):
    provider = ScriptedProvider([
        [tool_use("list_project", "t1")],
        [text("Annotating both genomes now."), tool_use("start_annotation", "t2", files=["demo_phage.fasta", "demo_rotated.fasta"])],
        [text("I started annotation of demo_phage and demo_rotated. I will report the results once they are complete.")],
    ])
    result = agent.run_turn(project, "Annotate these two genomes", provider)
    assert [c["name"] for c in result["tool_calls"]] == ["list_project", "start_annotation"]
    assert result["grounding_problems"] == []
    assert len([j for j in project.jobs() if j["kind"] == "annotation"]) == 2
    sent = provider.calls[-1]["messages"]
    assert sent[-1]["content"][0]["type"] == "tool_result"
    assert "Completed samples: none" in provider.calls[0]["system"]
    assert (project.root / "outbound.jsonl").is_file()


def test_agent_retries_then_flags_ungrounded_answer(project):
    provider = ScriptedProvider([[text("Annotation completed: PM_000777 is a capsid.")],
                                 [text("Annotation completed: PM_000777 is a capsid.")]])
    result = agent.run_turn(project, "What did you find?", provider)
    assert "Automatic grounding check failed" in json.dumps(provider.calls[1]["messages"][-1])
    assert result["grounding_problems"] and "Automatic check" in result["text"]


def test_agent_tool_error_is_returned_to_model_not_hidden(project):
    provider = ScriptedProvider([[tool_use("get_run_summary", "t1", sample="nope")],
                                 [text("There is no completed annotation for that sample yet.")]])
    result = agent.run_turn(project, "Summarise nope", provider)
    assert result["tool_calls"][0]["ok"] is False
    blocks = provider.calls[1]["messages"][-1]["content"]
    assert blocks[0]["is_error"] is True


def test_job_event_reports_once_after_all_jobs_end(compared):
    p, _, _ = compared
    event = agent.job_event(p)
    assert event and "COMPLETE" in event and "not a message from the researcher" in event
    assert agent.job_event(p) is None


def test_event_turns_do_not_count_as_approval(compared):
    p, _, _ = compared
    turns = p.user_turns()
    agent.run_turn(p, "[Automatic update] jobs ended", ScriptedProvider([[text("Both runs are COMPLETE.")]]), event=True)
    assert p.user_turns() == turns


def test_old_tool_results_are_shortened_for_the_model():
    messages = [{"role": "user", "content": [{"type": "tool_result", "tool_use_id": str(i), "content": "x" * 50}]} for i in range(20)]
    shortened = agent._for_model(messages)
    assert "shortened" in shortened[0]["content"][0]["content"]
    assert shortened[-1]["content"][0]["content"] == "x" * 50
    assert messages[0]["content"][0]["content"] == "x" * 50  # stored history untouched


def test_every_tool_has_a_valid_schema():
    for schema in tools.schemas():
        assert schema["input_schema"]["type"] == "object"
        assert set(schema["input_schema"]["required"]) <= set(schema["input_schema"]["properties"])
        assert schema["description"]
    names = {s["name"] for s in tools.schemas()}
    assert not names & {"shell", "run_code", "execute", "python"}
