"""Persistent metadata-only registry for locally installed evidence resources."""
from __future__ import annotations

import hashlib
import json
import os
import platform
import shutil
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
    PMFDB = "PMFDB"
    INPHARED_GENOMES = "INPHARED_GENOMES"
    REFSEQ = "REFSEQ"
    PHAGE_PROTEINS = "PHAGE_PROTEINS"
    CUSTOM = "CUSTOM"


class ResourceStatus(str, Enum):
    REGISTERED = "REGISTERED"
    READY = "READY"
    UNAVAILABLE = "UNAVAILABLE"
    INVALID = "INVALID"


REQUIRED_TOOLS = {
    "PFAM": "hmmscan", "VOGDB": "hmmscan", "SWISSPROT": "diamond",
    "PHROGS": "mmseqs", "PMFDB": "mmseqs", "INPHARED_GENOMES": "mash",
}


def default_registry_path() -> Path:
    override = os.environ.get("PHAGEMINE_REGISTRY_PATH")
    if override:
        return Path(override).expanduser()
    system = platform.system()
    if system == "Darwin":
        return Path.home() / "Library" / "Application Support" / "PhageMine" / "resources.json"
    if system == "Windows":
        return Path(os.environ.get("APPDATA", Path.home() / "AppData" / "Roaming")) / "PhageMine" / "resources.json"
    return Path(os.environ.get("XDG_CONFIG_HOME", Path.home() / ".config")) / "phagemine" / "resources.json"


@dataclass
class EvidenceResource:
    name: str
    resource_type: ResourceType | str
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
        data["resource_type"] = self.resource_type.value if isinstance(self.resource_type, ResourceType) else str(self.resource_type)
        data["status"] = self.status(check_checksum).value
        data["registry_checksum_policy"] = "checksum is optional and is not computed automatically"
        return data


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _tool_available(tool: str | None) -> bool:
    if not tool:
        return False
    candidate = Path(str(tool)).expanduser()
    return (candidate.is_file() and os.access(candidate, os.X_OK)) or bool(shutil.which(str(tool)))


def validate_resource(resource: dict[str, Any], *, check_checksum: bool = True) -> dict[str, Any]:
    """Validate a registered resource and its operational sidecars.

    This is deliberately conservative: a resource is READY only when the
    primary path, declared annotations/metadata, required executable(s), and
    format-specific index markers are present.
    """
    result = dict(resource)
    path = Path(str(resource.get("path", ""))).expanduser()
    problems: list[str] = []
    if not path.exists():
        result["status"] = ResourceStatus.UNAVAILABLE.value
        result["validation_errors"] = ["resource path does not exist"]
        return result
    if not os.access(path, os.R_OK):
        problems.append("resource path is not readable")
    if check_checksum and resource.get("checksum") and path.is_file() and _sha256(path) != resource["checksum"]:
        problems.append("resource checksum mismatch")
    provenance = resource.get("provenance") or {}
    kind = str(resource.get("resource_type", "")).upper()
    required_tool = REQUIRED_TOOLS.get(kind)
    declared_tools = list(resource.get("required_tools") or [])
    if required_tool and required_tool not in declared_tools:
        declared_tools.append(required_tool)
    for key in ("annotations_path", "metadata_path", "hmm_profiles_path"):
        value = provenance.get(key)
        if value and not Path(value).expanduser().is_file():
            problems.append(f"{key} does not exist")
    for tool in declared_tools:
        if not _tool_available(tool):
            problems.append(f"required executable unavailable: {tool}")
    if kind in {"PFAM", "VOGDB"} and path.is_file():
        # hmmsearch/hmmscan databases are operationally prepared when all
        # HMMER pressed companions exist; permit explicit small test fixtures.
        companions = [Path(str(path) + suffix) for suffix in (".h3f", ".h3i", ".h3m", ".h3p")]
        if not all(p.exists() for p in companions):
            problems.append("HMM database is marked prepared but pressed indexes are incomplete")
        if kind == "VOGDB" and not provenance.get("annotations_path"):
            problems.append("VOGDB annotation mapping is not registered")
    if kind == "PHROGS":
        if not Path(str(path) + ".dbtype").is_file():
            problems.append("MMseqs2 database .dbtype companion is missing")
        if not provenance.get("annotations_path"):
            problems.append("PHROGs annotation mapping is not registered")
    if kind == "SWISSPROT" and not provenance.get("metadata_path"):
        problems.append("Swiss-Prot metadata is not registered")
    if kind == "PMFDB":
        root = path if path.is_dir() else path.parent
        required = (
            "reference_phage_proteins.faa", "reference_metadata.tsv",
            "reference_qc.tsv", "reference_manifest.json",
            "mmseqs/target_db.dbtype", "mmseqs/mmseqs_index_manifest.json",
        )
        for name in required:
            if not (root / name).is_file():
                problems.append(f"PMFDB required file is missing: {name}")
        manifest_path = root / "reference_manifest.json"
        if manifest_path.is_file():
            try:
                manifest = json.loads(manifest_path.read_text())
                if str(manifest.get("schema_version")) != "1.0":
                    problems.append("PMFDB manifest schema_version is unsupported")
                for name, expected in (manifest.get("checksums") or {}).items():
                    candidate = root / name
                    if candidate.is_file() and _sha256(candidate) != expected:
                        problems.append(f"PMFDB checksum mismatch: {name}")
            except (OSError, ValueError, TypeError):
                problems.append("PMFDB reference_manifest.json is invalid")
    if kind == "INPHARED_GENOMES":
        root = path.parent if path.is_file() else path
        required = (
            "reference_phage_genomes.fna", "genome_metadata.tsv", "genome_qc.tsv",
            "genome_manifest.json", "inphared.msh",
        )
        for name in required:
            if not (root / name).is_file():
                problems.append(f"INPHARED required file is missing: {name}")
        manifest_path = root / "genome_manifest.json"
        if manifest_path.is_file():
            try:
                manifest = json.loads(manifest_path.read_text())
                if str(manifest.get("schema_version")) != "1.0":
                    problems.append("INPHARED manifest schema_version is unsupported")
                for name, expected in (manifest.get("checksums") or {}).items():
                    candidate = root / name
                    if candidate.is_file() and _sha256(candidate) != expected:
                        problems.append(f"INPHARED checksum mismatch: {name}")
            except (OSError, ValueError, TypeError):
                problems.append("INPHARED genome_manifest.json is invalid")
    result["status"] = ResourceStatus.INVALID.value if problems else ResourceStatus.READY.value
    result["validation_errors"] = problems
    return result


class EvidenceResourceManager:
    def __init__(self, registry_path: str | Path | None = None):
        self.registry_path = Path(registry_path).expanduser() if registry_path else default_registry_path()

    def _load(self) -> dict[str, EvidenceResource]:
        if not self.registry_path.exists():
            return {}
        payload = json.loads(self.registry_path.read_text())
        resources: dict[str, EvidenceResource] = {}
        for name, item in payload.get("resources", {}).items():
            raw_type = str(item["resource_type"]).upper()
            try:
                resource_type: ResourceType | str = ResourceType(raw_type)
            except ValueError:
                # Preserve forward compatibility: an entry written by a newer
                # PhageMine must not prevent older, known resources from loading.
                resource_type = raw_type
            resources[name] = EvidenceResource(
                name=name,
                resource_type=resource_type,
                **{key: value for key, value in item.items()
                   if key not in {"name", "resource_type", "status", "registry_checksum_policy"}},
            )
        return resources

    def _save(self, resources: dict[str, EvidenceResource]) -> None:
        self.registry_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {"format_version": "1", "resources": {name: resource.metadata() for name, resource in resources.items()}}
        encoded = json.dumps(payload, indent=2, sort_keys=True) + "\n"
        temporary = self.registry_path.with_name(self.registry_path.name + ".tmp")
        backup = self.registry_path.with_name(self.registry_path.name + ".bak")
        temporary.write_text(encoded)
        if self.registry_path.exists():
            shutil.copy2(self.registry_path, backup)
        os.replace(temporary, self.registry_path)

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

    def validate_all(self, check_checksum: bool = True) -> list[dict[str, Any]]:
        return [validate_resource(resource.metadata(check_checksum), check_checksum=check_checksum)
                for resource in self._load().values()]

    def validate(self, name: str, check_checksum: bool = True) -> dict[str, Any] | None:
        resource = self._load().get(name)
        return validate_resource(resource.metadata(check_checksum), check_checksum=check_checksum) if resource else None

    def get(self, name: str, check_checksum: bool = False) -> dict[str, Any] | None:
        resource = self._load().get(name)
        return resource.metadata(check_checksum) if resource else None

    def find(self, resource_type: ResourceType | str, check_checksum: bool = False) -> dict[str, Any] | None:
        target = resource_type if isinstance(resource_type, ResourceType) else ResourceType(str(resource_type).upper())
        # Legacy adapter lookup remains path-based for backwards compatibility;
        # strict operational selection is performed by preflight_profile and
        # validate_all before public execution.
        candidates = [resource.metadata(check_checksum) for resource in self._load().values()
                      if resource.resource_type == target and resource.status(check_checksum) == ResourceStatus.READY]
        return candidates[0] if candidates else None

    @staticmethod
    def checksum(path: str | Path) -> str:
        return _sha256(Path(path).expanduser())
