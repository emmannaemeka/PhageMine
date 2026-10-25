"""Local GUI deployment checks without opening a real browser."""
from pathlib import Path
import sys
import threading
from unittest.mock import MagicMock, patch

import pytest
from phagemine.gui import launcher


def test_browser_waits_for_health_before_opening():
    response = MagicMock()
    response.__enter__.return_value.status = 200
    opener = MagicMock()
    opener.open.side_effect = [OSError("not ready"), response]
    with patch.object(launcher.urllib.request, "build_opener", return_value=opener), patch.object(launcher.webbrowser, "open", return_value=True) as browser:
        launcher.open_when_ready("http://127.0.0.1:8507/", threading.Event())
    assert opener.open.call_count == 2
    browser.assert_called_once_with("http://127.0.0.1:8507/", new=2)


def test_stopped_server_does_not_open_browser():
    stopped = threading.Event(); stopped.set()
    with patch.object(launcher.webbrowser, "open") as browser:
        launcher.open_when_ready("http://127.0.0.1:8501/", stopped)
        browser.assert_not_called()


def test_browser_failure_prints_manual_link(capsys):
    response = MagicMock()
    response.__enter__.return_value.status = 200
    with patch.object(launcher.urllib.request, "build_opener") as factory, patch.object(launcher.webbrowser, "open", return_value=False):
        factory.return_value.open.return_value = response
        launcher.open_when_ready("http://127.0.0.1:8507/", threading.Event())
    assert "http://127.0.0.1:8507/" in capsys.readouterr().err


def test_busy_default_selects_another_port():
    server = MagicMock()
    server.__enter__.return_value = server
    server.bind.side_effect = [OSError("busy"), None]
    with patch.object(launcher.socket, "socket", return_value=server):
        assert launcher.available_port(None) == 8502
    assert server.bind.call_args.args == (("127.0.0.1", 8502),)


def test_busy_explicit_port_is_not_silently_changed():
    with patch.object(launcher.socket, "socket", side_effect=OSError("busy")):
        with pytest.raises(OSError, match="--port"):
            launcher.available_port(8507)


@pytest.mark.parametrize("args", [[], ["--no-browser"], ["--server.headless=true"]])
def test_browser_policy_and_local_server_arguments(args, monkeypatch):
    from streamlit.web import cli
    original = sys.argv
    commands = []
    monkeypatch.setattr(launcher, "available_port", lambda requested: 8507)
    monkeypatch.setattr(cli, "main", lambda: commands.append(list(sys.argv)) or 0)
    with patch.object(launcher.threading, "Thread") as worker:
        assert launcher.main(args) == 0
        assert worker.call_count == (0 if args else 1)
    assert sys.argv is original
    assert "--server.port=8507" in commands[0]
    assert "--server.address=127.0.0.1" in commands[0]
    assert "--server.headless=true" in commands[0]


def test_cli_gui_forwards_deployment_options(monkeypatch):
    from phagemine.cli import main
    calls = []
    monkeypatch.setattr(launcher, "main", lambda args: calls.append(args) or 0)
    assert main(["gui", "--port", "8510", "--no-browser"]) == 0
    assert calls == [["--port", "8510", "--no-browser"]]


def test_gui_is_required_by_default_installation():
    try:
        import tomllib
    except ImportError:
        import tomli as tomllib
    project = tomllib.loads((Path(__file__).parents[1] / "pyproject.toml").read_text())["project"]
    assert any(dep.startswith("streamlit") for dep in project["dependencies"])
    assert "gui" in project["optional-dependencies"]


def test_streamed_process_output_and_saved_log(tmp_path):
    from phagemine.gui.services.execution import execute
    lines = []
    result = execute([sys.executable, "-u", "-c", "import sys; print('starting', flush=True); print('finished', file=sys.stderr)", "--output", str(tmp_path)], on_output=lines.append)
    assert result.exit_code == 0
    assert "starting" in "".join(lines)
    assert "finished" in "".join(lines)
    assert (tmp_path / "gui_console.log").read_text() == result.stdout
    assert (tmp_path / "gui_metadata.json").is_file()
