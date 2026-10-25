"""Free model providers (Ollama, OpenAI-compatible) against fake HTTP servers."""
from __future__ import annotations

import json
import threading
from http.server import BaseHTTPRequestHandler, HTTPServer

import pytest

from phagemine.assistant import llm


class FakeServer:
    """Records requests; replies from a queue of canned responses."""

    def __init__(self, models=("qwen3:8b",)):
        self.requests, self.replies, self.models = [], [], list(models)
        owner = self

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass

            def _send(self, payload, code=200):
                body = json.dumps(payload).encode()
                self.send_response(code)
                self.send_header("Content-Type", "application/json")
                self.send_header("Content-Length", str(len(body)))
                self.end_headers()
                self.wfile.write(body)

            def do_GET(self):
                self._send({"models": [{"name": m} for m in owner.models]})

            def do_POST(self):
                length = int(self.headers["Content-Length"])
                owner.requests.append({"path": self.path, "headers": dict(self.headers),
                                       "body": json.loads(self.rfile.read(length))})
                self._send(owner.replies.pop(0))

        self.server = HTTPServer(("127.0.0.1", 0), Handler)
        self.url = f"http://127.0.0.1:{self.server.server_address[1]}"
        threading.Thread(target=self.server.serve_forever, daemon=True).start()

    def close(self):
        self.server.shutdown()


@pytest.fixture
def server():
    fake = FakeServer()
    yield fake
    fake.close()


TOOLS = [{"name": "query_proteins", "description": "Find proteins",
          "input_schema": {"type": "object", "properties": {"text": {"type": "string"}}, "required": []}}]

HISTORY = [
    {"role": "user", "content": "Find capsid genes"},
    {"role": "assistant", "content": [{"type": "text", "text": "Searching."},
                                      {"type": "tool_use", "id": "call_1", "name": "query_proteins", "input": {"text": "capsid"}}]},
    {"role": "user", "content": [{"type": "tool_result", "tool_use_id": "call_1", "content": "{\"total_matches\": 0}"}]},
]


def test_ollama_picks_best_installed_model_and_converts_tool_calls(server, monkeypatch):
    monkeypatch.delenv("PHAGEMINE_ASSISTANT_MODEL", raising=False)
    server.models = ["llama3.1:8b", "qwen3:14b", "nomic-embed-text"]
    provider = llm.OllamaProvider(base_url=server.url)
    assert provider.model == "qwen3:14b"
    server.replies.append({"message": {"role": "assistant", "content": "<think>plan</think>",
                                       "tool_calls": [{"function": {"name": "query_proteins", "arguments": {"text": "tail"}}}]},
                           "prompt_eval_count": 900, "eval_count": 20})
    result = provider.complete("SYSTEM", HISTORY, TOOLS)
    sent = server.requests[-1]["body"]
    assert server.requests[-1]["path"] == "/api/chat"
    assert sent["options"]["num_ctx"] >= 16384 and sent["stream"] is False
    assert sent["messages"][0] == {"role": "system", "content": "SYSTEM"}
    assert sent["messages"][2]["tool_calls"][0]["function"]["arguments"] == {"text": "capsid"}
    assert sent["messages"][3]["role"] == "tool" and sent["messages"][3]["name"] == "query_proteins"
    assert sent["tools"][0]["function"]["parameters"]["type"] == "object"
    assert result["stop_reason"] == "tool_use"
    assert [b["type"] for b in result["content"]] == ["tool_use"]  # <think> text removed
    assert result["content"][0]["input"] == {"text": "tail"} and result["content"][0]["id"].startswith("call_")


def test_ollama_reports_missing_model_and_server(server):
    server.models = ["nomic-embed-text"]
    with pytest.raises(llm.ProviderUnavailable, match="ollama pull qwen3:8b"):
        llm.OllamaProvider(base_url=server.url)
    with pytest.raises(llm.ProviderUnavailable, match="not installed"):
        llm.OllamaProvider(base_url=server.url, model="qwen3:32b")
    with pytest.raises(llm.ProviderUnavailable, match="not running"):
        llm.OllamaProvider(base_url="http://127.0.0.1:9")


def test_openai_compatible_round_trip(server):
    provider = llm.OpenAICompatibleProvider(base_url=server.url + "/v1", model="free-model", api_key="secret")
    server.replies.append({"choices": [{"message": {"role": "assistant", "content": None, "tool_calls": [
        {"id": "abc", "type": "function", "function": {"name": "query_proteins", "arguments": "{\"text\": \"portal\"}"}}]}}],
        "usage": {"prompt_tokens": 1000, "completion_tokens": 15}})
    result = provider.complete("SYSTEM", HISTORY, TOOLS)
    request = server.requests[-1]
    assert request["path"] == "/v1/chat/completions"
    assert request["headers"]["Authorization"] == "Bearer secret"
    assert request["body"]["messages"][2]["tool_calls"][0]["function"]["arguments"] == "{\"text\": \"capsid\"}"
    assert request["body"]["messages"][3] == {"role": "tool", "tool_call_id": "call_1", "name": "query_proteins",
                                              "content": "{\"total_matches\": 0}"}
    assert result["content"] == [{"type": "tool_use", "id": "abc", "name": "query_proteins", "input": {"text": "portal"}}]
    server.replies.append({"choices": [{"message": {"role": "assistant", "content": "No portal proteins were found."}}]})
    final = provider.complete("SYSTEM", HISTORY, TOOLS)
    assert final["stop_reason"] == "end_turn" and final["content"][0]["text"].startswith("No portal")


def test_malformed_tool_arguments_reach_the_tool_layer_as_errors(server):
    provider = llm.OpenAICompatibleProvider(base_url=server.url, model="m")
    server.replies.append({"choices": [{"message": {"tool_calls": [{"id": "x", "function": {"name": "query_proteins", "arguments": "{oops"}}]}}]})
    block = provider.complete("S", HISTORY, TOOLS)["content"][0]
    from phagemine.assistant import tools
    assert "_invalid_arguments" in block["input"]
    assert "unexpected" in tools.execute(tools.ToolContext(project=None), "query_proteins", block["input"])["error"]


def test_default_provider_selection(server, monkeypatch):
    for name in ("PHAGEMINE_ASSISTANT_PROVIDER", "ANTHROPIC_API_KEY", "PHAGEMINE_ASSISTANT_MODEL"):
        monkeypatch.delenv(name, raising=False)
    monkeypatch.setenv("OLLAMA_HOST", server.url)
    assert llm.default_provider().name == "ollama:qwen3:8b"
    monkeypatch.setenv("OLLAMA_HOST", "http://127.0.0.1:9")
    with pytest.raises(llm.ProviderUnavailable, match="Free options"):
        llm.default_provider()
    monkeypatch.setenv("PHAGEMINE_ASSISTANT_PROVIDER", "openai-compatible")
    monkeypatch.setenv("PHAGEMINE_ASSISTANT_BASE_URL", server.url)
    monkeypatch.setenv("PHAGEMINE_ASSISTANT_MODEL", "free-model")
    assert llm.default_provider().name == "openai-compatible:free-model"
    monkeypatch.setenv("PHAGEMINE_ASSISTANT_PROVIDER", "nonsense")
    with pytest.raises(llm.ProviderUnavailable, match="Unknown"):
        llm.default_provider()


def test_local_requests_bypass_proxies(server, monkeypatch):
    monkeypatch.setenv("HTTP_PROXY", "http://127.0.0.1:9")
    monkeypatch.setenv("http_proxy", "http://127.0.0.1:9")
    assert llm.OllamaProvider(base_url=server.url).model == "qwen3:8b"
