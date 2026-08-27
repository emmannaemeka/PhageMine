"""Small stderr-only pipeline progress reporter."""
from __future__ import annotations

import sys
import time
from typing import TextIO


class ProgressReporter:
    STAGES = ("input/genome validation", "gene prediction", "Pfam", "VOGDB", "Swiss-Prot", "PHROGs",
              "evidence integration", "candidate ranking/mining", "QC/report generation",
              "GenBank pre-submission package")

    def __init__(self, stream: TextIO | None = None, quiet: bool = False, no_progress: bool = False, callback=None):
        self.stream = stream or sys.stderr
        self.quiet = quiet
        self.no_progress = no_progress
        self.completed = 0
        self._started: tuple[str, float] | None = None
        self.callback = callback

    def start(self, stage: str) -> None:
        self._started = (stage, time.monotonic())
        if self.callback: self.callback("RUNNING", stage, None)
        if self.quiet:
            return
        self._write(f"[PhageMine {self.completed}/{len(self.STAGES)}] RUNNING: {stage}")

    def finish(self, summary: str | None = None) -> None:
        if self._started is None:
            return
        stage, started = self._started
        elapsed = time.monotonic() - started
        suffix = f"; {summary}" if summary else ""
        self.completed += 1
        if not self.quiet:
            self._write(f"[PhageMine {self.completed}/{len(self.STAGES)}] DONE: {stage} ({elapsed:.1f}s){suffix}")
        if self.callback: self.callback("COMPLETE", stage, {"elapsed_seconds": elapsed, "summary": summary})
        self._started = None

    def skip(self, summary: str | None = None) -> None:
        """Complete the current stage explicitly as unavailable/not configured."""
        if self._started is None:
            return

        stage, started = self._started
        elapsed = time.monotonic() - started
        suffix = f"; {summary}" if summary else ""

        self.completed += 1

        if not self.quiet:
            self._write(
                f"[PhageMine {self.completed}/{len(self.STAGES)}] "
                f"SKIPPED: {stage} ({elapsed:.1f}s){suffix}"
            )

        if self.callback:
            self.callback(
                "SKIPPED",
                stage,
                {"elapsed_seconds": elapsed, "summary": summary},
            )

        self._started = None

    def _write(self, message: str) -> None:
        self.stream.write(message + "\n")
        self.stream.flush()
