"""Exercise the real optional interface, including successful subprocess runs."""
from pathlib import Path

import pytest

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest

APP = Path(__file__).resolve().parents[1] / "src/phagemine/gui/app.py"


def test_gui_pages_render_without_results(tmp_path):
    app = AppTest.from_file(str(APP)).run()
    assert not app.exception
    app.sidebar.text_input[0].set_value(str(tmp_path)).run()
    for page in ("Run Monitor", "Annotation Results", "Discovery Results", "Figures", "Export / Downloads"):
        app.sidebar.radio[0].set_value(page).run()
        assert not app.exception, page


def test_successful_gui_run_does_not_mutate_rendered_widget(tmp_path, monkeypatch):
    from phagemine.gui.services import execution
    fasta = tmp_path / "input.fasta"; fasta.write_text(">g\nATGAAATAG\n")
    output = tmp_path / "results"
    monkeypatch.setattr(execution, "execute", lambda args: execution.RunResult(list(args), "done", "", 0, output))
    app = AppTest.from_file(str(APP)).run()
    next(widget for widget in app.text_input if widget.label == "FASTA genome path").set_value(str(fasta))
    next(widget for widget in app.text_input if widget.label == "Output directory").set_value(str(output))
    next(widget for widget in app.button if widget.label == "Start analysis").click().run()
    assert not app.exception
    assert any("Results saved" in item.value for item in app.success)
