"""Portable executable launch selection; never retry a failed analysis command."""
from __future__ import annotations

from functools import lru_cache
import platform
import signal
import subprocess


def _apple_silicon() -> bool:
    if platform.system() != "Darwin":
        return False
    if platform.machine() == "arm64":
        return True
    try:
        result = subprocess.run(
            ["/usr/sbin/sysctl", "-n", "hw.optional.arm64"],
            capture_output=True, text=True, timeout=5, check=False,
        )
        return result.returncode == 0 and result.stdout.strip() == "1"
    except (OSError, subprocess.SubprocessError):
        return False


@lru_cache(maxsize=64)
def mmseqs_prefix(executable: str) -> tuple[str, ...]:
    """Select ARM only after SIGILL and a successful ARM version probe.

    Selection is cached for this process. Pure Intel binaries are never marked
    compatible merely because they are executable. Linux/Windows are unchanged.
    """
    normal = (executable,)
    if not _apple_silicon():
        return normal
    try:
        result = subprocess.run(
            [executable, "version"], capture_output=True, text=True,
            timeout=10, check=False,
        )
        if result.returncode != -signal.SIGILL:
            return normal
        arm = ("/usr/bin/arch", "-arm64", executable)
        result = subprocess.run(
            [*arm, "version"], capture_output=True, text=True,
            timeout=10, check=False,
        )
        if result.returncode == 0 and (result.stdout or result.stderr).strip():
            return arm
    except (OSError, subprocess.SubprocessError):
        pass
    return normal


def run_mmseqs(command, **kwargs):
    """Launch version checks and work through the same selected architecture."""
    return subprocess.run([*mmseqs_prefix(str(command[0])), *command[1:]], **kwargs)
