"""Exercise the bundled interface, including successful subprocess runs."""
from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest

APP = Path(__file__).resolve().parents[1] / "src/phagemine/gui/app.py"


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
    app = AppTest.from_file(str(APP), default_timeout=15).run()
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
    app = AppTest.from_file(str(APP), default_timeout=15).run()
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
    app = AppTest.from_file(str(APP), default_timeout=15).run()
    next(w for w in app.text_input if w.label == "FASTA genome path").set_value(str(fasta))
    next(w for w in app.button if w.label == "Start analysis").click().run()
    assert not app.exception
    assert any("exited with code 2" in item.value for item in app.error)
    assert any("dependency could not start" in code.value for code in app.code)
    assert "last_completed_run" not in app.session_state


def test_assistant_search_and_scope_change(tmp_path):
    (tmp_path / 'annotation.tsv').write_text('protein_id\tproduct\tstart\tend\nPM_1\tmajor capsid protein\t1\t99\n')
    app = AppTest.from_file(str(APP), default_timeout=15).run()
    app.sidebar.text_input[0].set_value(str(tmp_path))
    app.sidebar.radio[0].set_value('Ask PhageMine').run()
    app.chat_input[0].set_value('Find the capsid gene').run()
    assert not app.exception
    assert any('Found 1 matching' in m.value for m in app.markdown)
    assert len(app.dataframe) == 1
    app.sidebar.text_input[0].set_value(str(tmp_path / 'other')).run()
    assert not app.exception
    assert len(app.dataframe) == 0


def test_assistant_annotation_handoff():
    app = AppTest.from_file(str(APP), default_timeout=15).run()
    app.sidebar.radio[0].set_value('Ask PhageMine').run()
    app.chat_input[0].set_value('Annotate my genome').run()
    next(b for b in app.button if b.label == 'Upload genome and start annotation').click().run()
    assert not app.exception
    assert app.sidebar.radio[0].value == 'Home / New Analysis'
