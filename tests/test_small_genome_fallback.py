from pathlib import Path
from types import SimpleNamespace

import pytest

from phagemine.reconciliation import ProdigalPredictor


@pytest.mark.parametrize("length, mode", [(5386, "meta"), (19999, "meta"),
                                          (20000, "single"), (20001, "single")])
def test_prodigal_mode_is_selected_by_supported_length(monkeypatch, tmp_path, length, mode):
    commands = []

    def fake_run(command, **_kwargs):
        commands.append(command)
        Path(command[command.index("-o") + 1]).write_text(
            "contig\tProdigal\tCDS\t1\t90\t.\t+\t0\tID=1_1\n"
        )
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr("phagemine.reconciliation.subprocess.run", fake_run)
    predictor = ProdigalPredictor("/test/prodigal")
    models = predictor.predict(tmp_path / "input.fa", "A" * length)
    assert predictor.parameters()["selected_mode"] == mode
    assert models[0].options["selected_mode"] == mode
    if mode == "meta":
        assert commands[0][-4:] == ["-p", "meta", "-g", "11"]
    else:
        assert "-p" not in commands[0] and "-g" not in commands[0]


@pytest.mark.parametrize("length", [5386, 25000])
def test_genuine_prodigal_failure_is_preserved(monkeypatch, tmp_path, length):
    commands = []

    def fake_run(command, **_kwargs):
        commands.append(command)
        return SimpleNamespace(returncode=7, stdout="", stderr="genuine tool failure")

    monkeypatch.setattr("phagemine.reconciliation.subprocess.run", fake_run)
    with pytest.raises(RuntimeError, match="exit 7.*genuine tool failure"):
        ProdigalPredictor("/test/prodigal").predict(tmp_path / "input.fa", "A" * length)
    assert len(commands) == 1
