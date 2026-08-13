"""Persistent metadata-only registry for locally installed evidence resources."""
from __future__ import annotations

import hashlib
import json
import os
import platform
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any


class ResourceType(str, Enum):
    PFAM = "PFAM"
    VOGDB = "VOGDB"
    PHROGS = "PHROGS"
    SWISSPROT = "SWISSPROT"
    REFSEQ = "REFSEQ"
    PHAGE_PROTEINS = "PHAGE_PROTEINS"
    CUSTOM = "CUSTOM"


class ResourceStatus(str, Enum):
    READY = "READY"
    UNAVAILABLE = "UNAVAILABLE"
    INVALID = "INVALID"


def default_registry_path() -> Path:
    system = platform.system()
    if system == "Darwin":
        return Path.home() / "Library" / "Application Support" / "PhageMine" / "resources.json"
    if system == "Windows":
        return Path(os.environ.get("APPDATA", Path.home() / "AppData" / "Roaming")) / "PhageMine" / "resources.json"
    return Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config")) / "phagemine" / "resources.json"


@dataclass
class EvidenceResource:
    name: str
    resource_type: ResourceType
    path: str
    version: str | None = None
    checksum: str | None = None
    date_registered: str = field(default_factory=lambda: datetime.now(timezone.utc).isoformat())
    preparation_status: str = "not_checked"
    required_tools: list[str] = field(default_factory=list)
    notes: str | None = None
    provenance: dict[str, Any] = field(default_factory=dict)
    checksum_algorithm: str | None = None

    def status(self, check_checksum: bool = False) -> ResourceStatus:
        resource_path = Path(self.path).expanduser()
        if not resource_path.exists():
            return ResourceStatus.UNAVAILABLE
        if not resource_path.is_file() and not resource_path.is_dir():
            return ResourceStatus.INVALID
        if not os.access(resource_path, os.R_OK):
            return ResourceStatus.INVALID
        if check_checksum and self.checksum:
            if resource_path.is_file() and _sha256(resource_path) != self.checksum:
                return ResourceStatus.INVALID
        return ResourceStatus.READY

    def metadata(self, check_checksum: bool = False) -> dict[str, Any]:
        data = asdict(self)
        data["resource_type"] = self.resource_type.value
        data["status"] = self.status(check_checksum).value
        data["registry_checksum_policy"] = "checksum is optional and is not computed automatically"
        return data


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class EvidenceResourceManager:
    def __init__(self, registry_path: str | Path | None = None):
        self.registry_path = Path(registry_path).expanduser() if registry_path else default_registry_path()

    def _load(self) -> dict[str, EvidenceResource]:
        if not self.registry_path.exists():
            return {}
        payload = json.loads(self.registry_path.read_text())
        return {name: EvidenceResource(name=name, resource_type=ResourceType(item["resource_type"]), **{key: value for key, value in item.items() if key not in {"name", "resource_type", "status", "registry_checksum_policy"}}) for name, item in payload.get("resources", {}).items()}

    def _save(self, resources: dict[str, EvidenceResource]) -> None:
        self.registry_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {"format_version": "1", "resources": {name: resource.metadata() for name, resource in resources.items()}}
        self.registry_path.write_text(json.dumps(payload, indent=2, sort_keys=True))

    def register(self, name: str, resource_type: ResourceType | str, path: str | Path, version: str | None = None, checksum: str | None = None, required_tools: list[str] | None = None, preparation_status: str = "not_checked", notes: str | None = None, provenance: dict[str, Any] | None = None) -> EvidenceResource:
        resource_type = resource_type if isinstance(resource_type, ResourceType) else ResourceType(str(resource_type).upper())
        resource = EvidenceResource(name=name, resource_type=resource_type, path=str(Path(path).expanduser()), version=version, checksum=checksum, required_tools=required_tools or [], preparation_status=preparation_status, notes=notes, provenance=provenance or {}, checksum_algorithm="sha256" if checksum else None)
        resources = self._load()
        resources[name] = resource
        self._save(resources)
        return resource

    def unregister(self, name: str) -> bool:
        resources = self._load()
        if name not in resources:
            return False
        del resources[name]
        self._save(resources)
        return True

    def list(self, check_checksum: bool = False) -> list[dict[str, Any]]:
        return [resource.metadata(check_checksum) for resource in self._load().values()]

    def get(self, name: str, check_checksum: bool = False) -> dict[str, Any] | None:
        resource = self._load().get(name)
        return resource.metadata(check_checksum) if resource else None

    def find(self, resource_type: ResourceType | str, check_checksum: bool = False) -> dict[str, Any] | None:
        target = resource_type if isinstance(resource_type, ResourceType) else ResourceType(str(resource_type).upper())
        for resource in self._load().values():
            if resource.resource_type == target:
                metadata = resource.metadata(check_checksum)
                if metadata["status"] == ResourceStatus.READY.value:
                    return metadata
        return None

    @staticmethod
    def checksum(path: str | Path) -> str:
        return _sha256(Path(path).expanduser())
