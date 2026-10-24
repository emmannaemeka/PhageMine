"""Architecture selection tests use fake tools, never biological inputs."""
import signal
import subprocess
from unittest.mock import patch

import pytest
from phagemine import executables as ex


def result(code=0, output="version-test\n"):
    return subprocess.CompletedProcess([], code, output, "")


@pytest.fixture(autouse=True)
def reset_selection():
    ex.mmseqs_prefix.cache_clear()
    yield
    ex.mmseqs_prefix.cache_clear()


@pytest.mark.parametrize("system", ["Linux", "Windows"])
def test_other_platforms_launch_unchanged(system):
    with patch.object(ex.platform, "system", return_value=system), patch.object(ex.subprocess, "run", return_value=result()) as run:
        ex.run_mmseqs(["/tools/custom-mmseqs", "version"], check=False)
        run.assert_called_once_with(["/tools/custom-mmseqs", "version"], check=False)


def test_rosetta_detection():
    with patch.object(ex.platform, "system", return_value="Darwin"), patch.object(ex.platform, "machine", return_value="x86_64"), patch.object(ex.subprocess, "run", return_value=result(output="1\n")) as run:
        assert ex._apple_silicon()
        assert run.call_args.args[0] == ["/usr/sbin/sysctl", "-n", "hw.optional.arm64"]


def test_successful_arm_probe_is_used_for_work_and_cached():
    with patch.object(ex, "_apple_silicon", return_value=True), patch.object(ex.subprocess, "run", side_effect=[result(-signal.SIGILL, ""), result(), result(), result()]) as run:
        ex.run_mmseqs(["/tools/mmseqs", "createdb", "input", "output"], check=True)
        ex.run_mmseqs(["/tools/mmseqs", "version"], check=False)
        assert run.call_count == 4
        assert run.call_args_list[1].args[0] == ["/usr/bin/arch", "-arm64", "/tools/mmseqs", "version"]
        assert run.call_args_list[2].args[0] == ["/usr/bin/arch", "-arm64", "/tools/mmseqs", "createdb", "input", "output"]
        assert run.call_args_list[2].kwargs == {"check": True}
        assert run.call_args_list[3].args[0] == ["/usr/bin/arch", "-arm64", "/tools/mmseqs", "version"]


@pytest.mark.parametrize("code", [0, 1, -11])
def test_non_sigill_never_triggers_arm_retry(code):
    with patch.object(ex, "_apple_silicon", return_value=True), patch.object(ex.subprocess, "run", return_value=result(code)) as run:
        assert ex.mmseqs_prefix("/tools/mmseqs") == ("/tools/mmseqs",)
        assert run.call_count == 1


def test_intel_only_binary_keeps_failure():
    with patch.object(ex, "_apple_silicon", return_value=True), patch.object(ex.subprocess, "run", side_effect=[result(-4, ""), result(1, "Bad CPU type")]):
        assert ex.mmseqs_prefix("/tools/mmseqs") == ("/tools/mmseqs",)


@pytest.mark.parametrize("error", [OSError("missing"), subprocess.TimeoutExpired("version", 10)])
def test_probe_errors_preserve_normal_invocation(error):
    with patch.object(ex, "_apple_silicon", return_value=True), patch.object(ex.subprocess, "run", side_effect=error):
        assert ex.mmseqs_prefix("/tools/mmseqs") == ("/tools/mmseqs",)


def test_work_failure_is_not_retried():
    with patch.object(ex, "_apple_silicon", return_value=False), patch.object(ex.subprocess, "run", return_value=result(-4, "")) as run:
        assert ex.run_mmseqs(["/tools/mmseqs", "search"]).returncode == -4
        assert run.call_count == 1


def test_doctor_uses_same_launch_selection(tmp_path):
    from phagemine.preflight import executable_status
    exe = tmp_path / "mmseqs"
    exe.write_text("fake")
    exe.chmod(0o755)
    with patch.object(ex, "_apple_silicon", return_value=True), patch.object(ex.subprocess, "run", side_effect=[result(-4, ""), result(), result()]):
        state = executable_status("mmseqs", str(exe))
    assert state["status"] == "READY"
    assert state["launch_command"] == ["/usr/bin/arch", "-arm64", str(exe)]
