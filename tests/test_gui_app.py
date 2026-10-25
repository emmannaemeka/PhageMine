"""Exercise the bundled interface, including successful subprocess runs."""
from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest

APP = Path(__file__).resolve().parents[1] / "src/phagemine/gui/app.py"


@pytest.fixture(autouse=True)
def isolated_assistant(tmp_path, monkeypatch):
    # The assistant page is the default; keep its projects out of the user's home.
    monkeypatch.setenv("PHAGEMINE_PROJECTS", str(tmp_path / "assistant-projects"))
    monkeypatch.delenv("ANTHROPIC_API_KEY", raising=False)


def home_page():
    app = AppTest.from_file(str(APP), default_timeout=15).run()
    app.sidebar.radio[0].set_value("Home / New Analysis").run()
    return app


def test_gui_pages_render_without_results(tmp_path):
    app = AppTest.from_file(str(APP), default_timeout=15).run()
    assert not app.exception
    app.sidebar.text_input[0].set_value(str(tmp_path)).run()
    for page in ("Run Monitor", "Annotation Results", "Discovery Results", "Figures", "Export / Downloads"):
        app.sidebar.radio[0].set_value(page).run()
        assert not app.exception, page


def test_successful_gui_run_does_not_mutate_rendered_widget(tmp_path, monkeypatch):
    from phagemine.gui.services import execution
    fasta = tmp_path / "input.fasta"; fasta.write_text(">g\nATGAAATAG\n")
    output = tmp_path / "results"
    monkeypatch.setattr(execution, "execute", lambda args, **kwargs: execution.RunResult(list(args), "done", "", 0, output))
    app = home_page()
    next(widget for widget in app.text_input if widget.label == "FASTA genome path").set_value(str(fasta))
    next(widget for widget in app.text_input if widget.label == "Output directory").set_value(str(output))
    next(widget for widget in app.button if widget.label == "Start analysis").click().run()
    assert not app.exception
    assert any("Results saved" in item.value for item in app.success)


def test_completed_run_opens_selected_results(tmp_path, monkeypatch):
    from phagemine.gui.services import execution
    fasta = tmp_path / "input.fasta"; fasta.write_text(">g\nATG\n")
    output = tmp_path / "completed"
    monkeypatch.setattr(execution, "execute", lambda args, **kwargs: execution.RunResult(list(args), "done", "", 0, output))
    app = home_page()
    next(w for w in app.text_input if w.label == "FASTA genome path").set_value(str(fasta))
    next(w for w in app.text_input if w.label == "Output directory").set_value(str(output))
    next(w for w in app.button if w.label == "Start analysis").click().run()
    next(w for w in app.button if w.label == "View results").click().run()
    assert not app.exception
    assert app.sidebar.text_input[0].value == str(output)
    assert app.sidebar.radio[0].value == "Annotation Results"


def test_prepared_download_survives_rerun(tmp_path):
    (tmp_path / "report.txt").write_text("test report")
    app = AppTest.from_file(str(APP), default_timeout=15).run()
    app.sidebar.text_input[0].set_value(str(tmp_path))
    app.sidebar.radio[0].set_value("Export / Downloads").run()
    next(w for w in app.button if w.label == "Prepare result ZIP").click().run()
    assert not app.exception
    assert app.session_state["result_zip"][1].startswith(b"PK")
    app.run()
    assert not app.exception
    assert app.session_state["result_zip"][0] == str(tmp_path.resolve())


def test_environment_guidance_renders_without_running_external_tools(monkeypatch):
    from phagemine import database_health, preflight
    monkeypatch.setattr(database_health, "health", lambda **kwargs: {"missing_annotation_resources": [], "resources": []})
    monkeypatch.setattr(preflight, "doctor", lambda **kwargs: {"executables": [], "resources": [], "recommendations": [{"commands": ["phagemine doctor"]}]})
    app = AppTest.from_file(str(APP), default_timeout=15).run()
    app.sidebar.radio[0].set_value("Environment / Database Status").run()
    assert not app.exception
    assert any("phagemine doctor" in code.value for code in app.code)


def test_failed_run_displays_live_diagnostic(tmp_path, monkeypatch):
    from phagemine.gui.services import execution
    fasta = tmp_path / "input.fasta"; fasta.write_text(">g\nATG\n")
    def fail(args, on_output):
        on_output("A dependency could not start")
        return execution.RunResult(list(args), "A dependency could not start", "", 2, tmp_path / "failed")
    monkeypatch.setattr(execution, "execute", fail)
    app = home_page()
    next(w for w in app.text_input if w.label == "FASTA genome path").set_value(str(fasta))
    next(w for w in app.button if w.label == "Start analysis").click().run()
    assert not app.exception
    assert any("exited with code 2" in item.value for item in app.error)
    assert any("dependency could not start" in code.value for code in app.code)
    assert "last_completed_run" not in app.session_state


def test_assistant_is_the_default_page_and_explains_missing_key():
    app = AppTest.from_file(str(APP), default_timeout=15).run()
    assert not app.exception
    assert app.sidebar.radio[0].value == "Ask PhageMine"
    assert any("Describe what you want to investigate" in item.value for item in app.markdown)
    assert any("ANTHROPIC_API_KEY" in item.value for item in app.error)


def test_assistant_project_creation_and_transcript_render(tmp_path):
    app = AppTest.from_file(str(APP), default_timeout=15).run()
    next(w for w in app.sidebar.text_input if w.label == "New project name").set_value("Phage pair")
    app.sidebar.button[0].click().run()
    assert not app.exception
    from phagemine.assistant.store import Project
    project = Project.list()[0]
    assert project.name == "Phage pair"
    project.append_transcript({"role": "user", "text": "Annotate my genome"})
    project.append_transcript({"role": "assistant", "text": "Started annotation.", "tool_calls": [{"name": "start_annotation", "input": {}, "ok": True}]})
    app.run()
    assert not app.exception
    assert any("Started annotation." in item.value for item in app.markdown)


class ReadingModel:
    """Stand-in model that decides from tool results, to drive the page end to end."""

    name = "test-reading-model"

    def complete(self, system, messages, tools):
        last = messages[-1]["content"]
        if isinstance(last, str):
            if "Automatic update" in last:
                return self._use("get_shared_genes", {})
            if "unknown" in last.lower():
                return self._use("query_proteins", {"unknown_function_only": True})
            return self._use("start_annotation", {"files": ["demo_phage.fasta", "demo_copy.fasta"]})
        import json
        result = json.loads(last[0]["content"])
        if "started" in result:
            ids = [item["job_id"] for item in result["started"]]
            return self._use("start_comparison", {"annotation_job_ids": ids, "allow_exact_only": True})
        if "waiting_for" in result:
            return self._say("Annotation and comparison are running. I will report the results once they are complete.")
        if "proteins_with_homolog_in_other" in result:
            counts = result["proteins_with_homolog_in_other"]
            return self._say("The comparison is COMPLETE. Proteins with a homolog in the other genome: "
                             + ", ".join(f"{k}: {v}" for k, v in counts.items()) + ".")
        if "proteins" in result:
            return self._say("Unknown-function proteins: " + ", ".join(p["protein_id"] for p in result["proteins"]) + ".")
        return self._say("Done.")

    def _use(self, name, arguments):
        return {"content": [{"type": "tool_use", "id": f"t{id(arguments)}", "name": name, "input": arguments}], "stop_reason": "tool_use", "usage": {}}

    def _say(self, words):
        return {"content": [{"type": "text", "text": words}], "stop_reason": "end_turn", "usage": {}}


def test_conversational_analysis_end_to_end(monkeypatch):
    import time
    from phagemine.assistant import tools
    from phagemine.assistant.store import Project
    monkeypatch.setenv("PHAGEMINE_ASSISTANT_DEMO", "1")
    monkeypatch.setattr(tools.shutil, "which", lambda name: None)
    demo = Path(__file__).resolve().parents[1] / "examples" / "demo_phage.fasta"
    project = Project.create("end to end")
    project.add_input("demo_phage.fasta", demo.read_bytes())
    project.add_input("demo_copy.fasta", demo.read_bytes().replace(b">demo_phage", b">demo_copy"))
    app = AppTest.from_file(str(APP), default_timeout=30).run()
    app.session_state["assistant_provider"] = ReadingModel()
    app.run()
    app.chat_input[0].set_value("Here are two phage genomes. Annotate them, compare them, and show me which genes they share.").run()
    assert not app.exception
    assert len(project.jobs()) == 3
    deadline = time.time() + 120
    while time.time() < deadline and not any("is COMPLETE" in m.value for m in app.markdown):
        time.sleep(1)
        app.run()
        assert not app.exception
    assert any("demo_phage: 5" in m.value for m in app.markdown), [m.value for m in app.markdown]
    app.chat_input[0].set_value("Which of these genes have unknown functions?").run()
    assert any("Unknown-function proteins: PM_" in m.value for m in app.markdown)
    transcript = Project(project.root).transcript()
    assert [t["role"] for t in transcript] == ["user", "assistant", "event", "assistant", "user", "assistant"]
    assert all(not t.get("grounding_problems") for t in transcript if t["role"] == "assistant")
