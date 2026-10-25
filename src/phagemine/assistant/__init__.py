"""Conversational assistant over PhageMine's tested workflows.

The language model plans and explains. PhageMine's existing CLI performs all
scientific computation, and deterministic code decides whether a job finished
and whether an answer is grounded in recorded tool outputs.
"""
from __future__ import annotations

ASSISTANT_VERSION = "0.1.0"

__all__ = ["ASSISTANT_VERSION"]
