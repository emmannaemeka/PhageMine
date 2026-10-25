"""Model providers behind one small interface.

``complete(system, messages, tools)`` returns ``{"content": [blocks], "stop_reason": str,
"usage": dict}`` where blocks are plain dicts in Anthropic Messages format.
A local or other provider only has to return the same shape.
"""
from __future__ import annotations

import os
from typing import Protocol

DEFAULT_MODEL = "claude-sonnet-5-5"


class Provider(Protocol):
    name: str

    def complete(self, system: str, messages: list[dict], tools: list[dict]) -> dict: ...


class ProviderUnavailable(RuntimeError):
    pass


class AnthropicProvider:
    """Claude through the Anthropic API. Needs ``pip install phagemine[assistant]`` and ANTHROPIC_API_KEY."""

    def __init__(self, model: str | None = None, api_key: str | None = None, max_tokens: int = 4096):
        try:
            import anthropic
        except ImportError as exc:  # pragma: no cover - depends on optional install
            raise ProviderUnavailable("The assistant needs the Anthropic SDK: pip install 'phagemine[assistant]'") from exc
        key = api_key or os.environ.get("ANTHROPIC_API_KEY")
        if not key:
            raise ProviderUnavailable("Set the ANTHROPIC_API_KEY environment variable to use the assistant")
        self.client = anthropic.Anthropic(api_key=key)
        self.model = model or os.environ.get("PHAGEMINE_ASSISTANT_MODEL") or DEFAULT_MODEL
        self.max_tokens = max_tokens
        self.name = f"anthropic:{self.model}"

    def complete(self, system: str, messages: list[dict], tools: list[dict]) -> dict:
        response = self.client.messages.create(
            model=self.model, max_tokens=self.max_tokens,
            system=[{"type": "text", "text": system, "cache_control": {"type": "ephemeral"}}],
            tools=tools, messages=messages)
        content = []
        for block in response.content:
            if block.type == "text":
                content.append({"type": "text", "text": block.text})
            elif block.type == "tool_use":
                content.append({"type": "tool_use", "id": block.id, "name": block.name, "input": block.input})
        usage = getattr(response, "usage", None)
        return {"content": content, "stop_reason": response.stop_reason,
                "usage": {"input_tokens": getattr(usage, "input_tokens", None),
                          "output_tokens": getattr(usage, "output_tokens", None)}}


class ScriptedProvider:
    """Replays fixed responses; used by tests to exercise the loop without a model."""

    name = "scripted"

    def __init__(self, responses: list[list[dict]]):
        self.responses = list(responses)
        self.calls: list[dict] = []

    def complete(self, system: str, messages: list[dict], tools: list[dict]) -> dict:
        self.calls.append({"system": system, "messages": [dict(m) for m in messages], "tools": tools})
        if not self.responses:
            return {"content": [{"type": "text", "text": "(no scripted response left)"}], "stop_reason": "end_turn", "usage": {}}
        content = self.responses.pop(0)
        stop = "tool_use" if any(b["type"] == "tool_use" for b in content) else "end_turn"
        return {"content": content, "stop_reason": stop, "usage": {}}


def default_provider() -> Provider:
    return AnthropicProvider()
