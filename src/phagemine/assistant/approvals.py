"""Researcher approvals that the model cannot grant.

Some actions need the researcher's explicit decision: applying NCBI corrections,
or comparing genomes without MMseqs2. Tools may only *request* an approval. The
chat page shows each pending request with Approve / Decline buttons, and only
those buttons call ``approve`` or ``decline``. An approval is bound to a hash of
the exact payload shown, so changing the payload after approval voids it.
"""
from __future__ import annotations

import hashlib
import json
import re
import uuid

from .store import Project, now

_ID = re.compile(r"approval-[0-9a-f]{8}")


def payload_hash(payload) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()


def _path(approval_id: str) -> str:
    if not _ID.fullmatch(approval_id or ""):
        raise ValueError("Unknown approval ID")
    return f"approvals/{approval_id}.json"


def request(project: Project, kind: str, payload, title: str) -> dict:
    record = {"id": f"approval-{uuid.uuid4().hex[:8]}", "kind": kind, "title": title, "payload": payload,
              "payload_hash": payload_hash(payload), "status": "PENDING", "created": now(), "used": False}
    project.write_state(_path(record["id"]), record)
    return record


def get(project: Project, approval_id: str) -> dict:
    record = project.read_state(_path(approval_id), None)
    if not record:
        raise ValueError("Unknown approval ID")
    return record


def _decide(project: Project, approval_id: str, status: str) -> dict:
    record = get(project, approval_id)
    if record["status"] != "PENDING":
        return record
    record.update(status=status, decided=now())
    project.write_state(_path(approval_id), record)
    return record


def approve(project: Project, approval_id: str) -> dict:
    """Called only from the researcher's Approve button, never from a tool."""
    return _decide(project, approval_id, "APPROVED")


def decline(project: Project, approval_id: str) -> dict:
    return _decide(project, approval_id, "DECLINED")


def pending(project: Project) -> list[dict]:
    folder = project.root / "approvals"
    if not folder.is_dir():
        return []
    records = [json.loads(p.read_text()) for p in folder.glob("approval-*.json")]
    return sorted((r for r in records if r["status"] == "PENDING"), key=lambda r: r["created"])


def check(project: Project, approval_id: str, kind: str, payload) -> None:
    """Raise unless this approval was granted for exactly this payload and not yet used."""
    record = get(project, approval_id)
    if record["kind"] != kind:
        raise ValueError("Approval is for a different action")
    if record["status"] == "DECLINED":
        raise ValueError("The researcher declined this request")
    if record["status"] != "APPROVED":
        raise ValueError("Not approved yet: the researcher must press Approve on the request shown in the chat page")
    if record["payload_hash"] != payload_hash(payload):
        raise ValueError("The request changed after it was approved; ask for approval again")
    if record.get("used"):
        raise ValueError("This approval was already used")


def consume(project: Project, approval_id: str) -> None:
    record = get(project, approval_id)
    record.update(used=True, used_at=now())
    project.write_state(_path(approval_id), record)


def unused_approved(project: Project, kind: str) -> dict | None:
    folder = project.root / "approvals"
    if not folder.is_dir():
        return None
    for path in sorted(folder.glob("approval-*.json")):
        record = json.loads(path.read_text())
        if record["kind"] == kind and record["status"] == "APPROVED" and not record.get("used"):
            return record
    return None
