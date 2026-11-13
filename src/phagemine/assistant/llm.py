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


import json as _json
import re as _re
import urllib.error as _urlerror
import urllib.request as _urlrequest
import uuid as _uuid
from urllib.parse import urlparse as _urlparse

OLLAMA_URL = "http://localhost:11434"
# Preferred local models, best first. Qwen3 models are trained for tool calling.
OLLAMA_PREFERRED = ("qwen3:32b", "qwen3:30b-a3b", "qwen3:14b", "qwen3:8b", "qwen3:4b",
                    "qwen2.5:14b", "qwen2.5:7b", "llama3.1:8b", "mistral-nemo")
_THINK = _re.compile(r"<think>.*?</think>", _re.S)


def _post(url: str, payload: dict, headers: dict | None = None, timeout: float = 900) -> dict:
    host = (_urlparse(url).hostname or "").lower()
    handlers = [_urlrequest.ProxyHandler({})] if host in {"localhost", "127.0.0.1", "::1"} else []
    opener = _urlrequest.build_opener(*handlers)
    request = _urlrequest.Request(url, data=_json.dumps(payload).encode(), method="POST",
                                  headers={"Content-Type": "application/json", **(headers or {})})
    try:
        with opener.open(request, timeout=timeout) as response:
            return _json.loads(response.read())
    except _urlerror.HTTPError as exc:
        detail = exc.read().decode(errors="replace")[:500]
        raise RuntimeError(f"Model server returned HTTP {exc.code}: {detail}") from exc
    except _urlerror.URLError as exc:
        raise ProviderUnavailable(f"Cannot reach the model server at {url}: {exc.reason}") from exc


def _arguments(raw) -> dict:
    if isinstance(raw, dict):
        return raw
    try:
        value = _json.loads(raw or "{}")
        return value if isinstance(value, dict) else {"_invalid_arguments": raw}
    except ValueError:
        # Passed through so the tool layer rejects it and the model sees why.
        return {"_invalid_arguments": raw}


def _chat_messages(system: str, messages: list[dict], *, arguments_as_text: bool) -> list[dict]:
    """Anthropic-format history -> OpenAI/Ollama chat messages."""
    names = {}
    out = [{"role": "system", "content": system}]
    for message in messages:
        content = message["content"]
        if isinstance(content, str):
            out.append({"role": message["role"], "content": content})
            continue
        if message["role"] == "assistant":
            text = "\n".join(b["text"] for b in content if b.get("type") == "text")
            calls = []
            for block in content:
                if block.get("type") == "tool_use":
                    names[block["id"]] = block["name"]
                    arguments = _json.dumps(block.get("input") or {}) if arguments_as_text else (block.get("input") or {})
                    calls.append({"id": block["id"], "type": "function",
                                  "function": {"name": block["name"], "arguments": arguments}})
            entry = {"role": "assistant", "content": text or ("" if arguments_as_text else "")}
            if calls:
                entry["tool_calls"] = calls
            out.append(entry)
        else:
            for block in content:
                if block.get("type") == "tool_result":
                    out.append({"role": "tool", "tool_call_id": block["tool_use_id"],
                                "name": names.get(block["tool_use_id"], ""), "content": str(block.get("content", ""))})
                elif block.get("type") == "text":
                    out.append({"role": "user", "content": block["text"]})
    return out


def _chat_tools(tools: list[dict]) -> list[dict]:
    return [{"type": "function", "function": {"name": t["name"], "description": t["description"],
                                              "parameters": t["input_schema"]}} for t in tools]


def _blocks(message: dict) -> list[dict]:
    blocks = []
    text = _THINK.sub("", message.get("content") or "").strip()
    if text:
        blocks.append({"type": "text", "text": text})
    for call in message.get("tool_calls") or []:
        function = call.get("function", {})
        blocks.append({"type": "tool_use", "id": call.get("id") or f"call_{_uuid.uuid4().hex[:12]}",
                       "name": function.get("name", ""), "input": _arguments(function.get("arguments"))})
    return blocks or [{"type": "text", "text": ""}]


class OllamaProvider:
    """A free local model through Ollama (https://ollama.com). Nothing leaves the computer."""

    def __init__(self, model: str | None = None, base_url: str | None = None, context: int | None = None):
        self.base_url = (base_url or os.environ.get("OLLAMA_HOST") or OLLAMA_URL).rstrip("/")
        if not self.base_url.startswith("http"):
            self.base_url = "http://" + self.base_url
        installed = self.installed_models(self.base_url)
        wanted = model or os.environ.get("PHAGEMINE_ASSISTANT_MODEL")
        if wanted:
            if wanted not in installed and f"{wanted}:latest" not in installed:
                raise ProviderUnavailable(f"Ollama model '{wanted}' is not installed. Run: ollama pull {wanted}")
            self.model = wanted
        else:
            self.model = next((m for m in OLLAMA_PREFERRED if m in installed), None)
            if not self.model:
                raise ProviderUnavailable("Ollama is running but has no suitable model. Run: ollama pull qwen3:8b")
        # Ollama's default context is too small for the tool list; ask for enough room.
        self.context = int(context or os.environ.get("PHAGEMINE_ASSISTANT_CONTEXT") or 16384)
        self.name = f"ollama:{self.model}"

    @staticmethod
    def installed_models(base_url: str = OLLAMA_URL, timeout: float = 2) -> list[str]:
        opener = _urlrequest.build_opener(_urlrequest.ProxyHandler({}))
        try:
            with opener.open(base_url.rstrip("/") + "/api/tags", timeout=timeout) as response:
                return [m["name"] for m in _json.loads(response.read()).get("models", [])]
        except (OSError, ValueError, _urlerror.URLError) as exc:
            raise ProviderUnavailable(f"Ollama is not running at {base_url}. Install it from https://ollama.com and start it.") from exc

    def complete(self, system: str, messages: list[dict], tools: list[dict]) -> dict:
        payload = {"model": self.model, "stream": False, "tools": _chat_tools(tools),
                   "messages": _chat_messages(system, messages, arguments_as_text=False),
                   "options": {"num_ctx": self.context, "temperature": 0.2}}
        data = _post(self.base_url + "/api/chat", payload)
        blocks = _blocks(data.get("message") or {})
        stop = "tool_use" if any(b["type"] == "tool_use" for b in blocks) else "end_turn"
        return {"content": blocks, "stop_reason": stop,
                "usage": {"input_tokens": data.get("prompt_eval_count"), "output_tokens": data.get("eval_count")}}


class OpenAICompatibleProvider:
    """Any server speaking the OpenAI chat-completions format: free cloud tiers
    (for example Google Gemini, Groq, OpenRouter) or a local server such as LM Studio."""

    def __init__(self, base_url: str | None = None, model: str | None = None, api_key: str | None = None):
        self.base_url = (base_url or os.environ.get("PHAGEMINE_ASSISTANT_BASE_URL") or "").rstrip("/")
        self.model = model or os.environ.get("PHAGEMINE_ASSISTANT_MODEL")
        self.api_key = api_key or os.environ.get("PHAGEMINE_ASSISTANT_API_KEY")
        if not self.base_url or not self.model:
            raise ProviderUnavailable("Set PHAGEMINE_ASSISTANT_BASE_URL and PHAGEMINE_ASSISTANT_MODEL "
                                      "(and PHAGEMINE_ASSISTANT_API_KEY if the service needs one)")
        self.name = f"openai-compatible:{self.model}"

    def complete(self, system: str, messages: list[dict], tools: list[dict]) -> dict:
        payload = {"model": self.model, "tools": _chat_tools(tools), "temperature": 0.2,
                   "messages": _chat_messages(system, messages, arguments_as_text=True)}
        headers = {"Authorization": f"Bearer {self.api_key}"} if self.api_key else {}
        data = _post(self.base_url + "/chat/completions", payload, headers)
        choice = (data.get("choices") or [{}])[0]
        blocks = _blocks(choice.get("message") or {})
        stop = "tool_use" if any(b["type"] == "tool_use" for b in blocks) else "end_turn"
        usage = data.get("usage") or {}
        return {"content": blocks, "stop_reason": stop,
                "usage": {"input_tokens": usage.get("prompt_tokens"), "output_tokens": usage.get("completion_tokens")}}


PROVIDERS = {"anthropic": AnthropicProvider, "ollama": OllamaProvider, "openai-compatible": OpenAICompatibleProvider}

SETUP_HELP = (
    "No language model is configured. Free options:\n"
    "1. Local and private: install Ollama (https://ollama.com), run `ollama pull qwen3:8b`, keep Ollama open, restart PhageMine.\n"
    "2. Free cloud tier: set PHAGEMINE_ASSISTANT_PROVIDER=openai-compatible with PHAGEMINE_ASSISTANT_BASE_URL, "
    "PHAGEMINE_ASSISTANT_MODEL and PHAGEMINE_ASSISTANT_API_KEY.\n"
    "Paid: set ANTHROPIC_API_KEY. See docs/assistant.md.")


def default_provider() -> Provider:
    """PHAGEMINE_ASSISTANT_PROVIDER if set; otherwise Anthropic when a key exists,
    then a running Ollama. Raises ProviderUnavailable with setup help."""
    choice = os.environ.get("PHAGEMINE_ASSISTANT_PROVIDER", "").strip().lower()
    if choice:
        if choice not in PROVIDERS:
            raise ProviderUnavailable(f"Unknown PHAGEMINE_ASSISTANT_PROVIDER '{choice}'. Use one of: {', '.join(PROVIDERS)}")
        return PROVIDERS[choice]()
    if os.environ.get("ANTHROPIC_API_KEY"):
        return AnthropicProvider()
    try:
        return OllamaProvider()
    except ProviderUnavailable as exc:
        if "no suitable model" in str(exc):
            raise
        raise ProviderUnavailable(SETUP_HELP) from exc
