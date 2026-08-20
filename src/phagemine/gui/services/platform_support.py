"""Read-only platform backend diagnostics for the desktop readiness screen."""
from __future__ import annotations

import platform
import shutil
import subprocess


def windows_wsl_status() -> dict:
    if platform.system() != "Windows":
        return {"applicable": False, "state": "NOT APPLICABLE", "detail": "WSL2 is only used on Windows."}
    executable = shutil.which("wsl.exe")
    if not executable:
        return {"applicable": True, "state": "NOT CONFIGURED",
                "detail": "WSL2 is not installed. Windows full scientific analysis is not ready."}
    try:
        result = subprocess.run([executable, "--status"], capture_output=True, text=True,
                                timeout=15, check=False, shell=False)
    except (OSError, subprocess.SubprocessError) as exc:
        return {"applicable": True, "state": "ERROR", "detail": str(exc)}
    text = (result.stdout or result.stderr).strip()
    state = "AVAILABLE" if result.returncode == 0 else "NOT READY"
    return {"applicable": True, "state": state,
            "detail": text[:1000] or f"wsl.exe exited with code {result.returncode}"}
