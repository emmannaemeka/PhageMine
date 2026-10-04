"""Stream self-contained, checksummed result archives and verify them offline."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path, PurePosixPath
import tempfile
import zipfile

from . import __version__
from .artifact_ops import checksum, completed_run


def bundle(run: str | Path, output: str | Path) -> dict:
    root = completed_run(run)
    output = Path(output).expanduser().resolve()
    if output == root or root in output.parents:
        raise ValueError("Write the archive outside its source directory")
    if output.exists(): raise ValueError("Archive already exists")
    files = sorted(path for path in root.rglob("*") if path.is_file())
    if any(path.is_symlink() for path in root.rglob("*")):
        raise ValueError("Symbolic links are not allowed in reproducibility bundles")
    if any(path.relative_to(root).as_posix() == "reproducibility_manifest.json" for path in files):
        files = [path for path in files if path.name != "reproducibility_manifest.json"]
    hashes = {path.relative_to(root).as_posix(): checksum(path) for path in files}
    manifest = {"schema_version": 1, "format": "PHAGEMINE_REPRODUCIBILITY_BUNDLE",
                "exporter_version": __version__, "hash_algorithm": "SHA256", "files": hashes,
                "includes_third_party_databases": False,
                "database_note": "Run and evidence manifests retain database identity; database bytes are not copied"}
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = tempfile.NamedTemporaryFile(prefix=".phagemine-bundle-", dir=output.parent, delete=False)
    temporary.close(); path = Path(temporary.name)
    try:
        with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED, allowZip64=True) as archive:
            for source in files:
                archive.write(source, source.relative_to(root).as_posix())
                if checksum(source) != hashes[source.relative_to(root).as_posix()]:
                    raise ValueError("Source artifacts changed during export; retry on a completed run")
            archive.writestr("reproducibility_manifest.json", json.dumps(manifest, indent=2, sort_keys=True) + "\n")
        verification = verify(path)
        if output.exists(): raise ValueError("Archive destination appeared during export")
        path.rename(output)
    finally:
        path.unlink(missing_ok=True)
    return {"status": "EXPORTED", "file_count": len(hashes), "archive_sha256": checksum(output), "output": str(output), "verification": verification}


def verify(path: str | Path) -> dict:
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        if len(names) != len(set(names)):
            raise ValueError("Duplicate archive entries")
        for name in names:
            part = PurePosixPath(name)
            if part.is_absolute() or ".." in part.parts or "\\" in name:
                raise ValueError("Unsafe archive entry")
        if "reproducibility_manifest.json" not in names:
            raise ValueError("Reproducibility manifest is missing")
        manifest = json.loads(archive.read("reproducibility_manifest.json"))
        if manifest.get("schema_version") != 1 or not isinstance(manifest.get("files"), dict):
            raise ValueError("Unsupported bundle manifest")
        expected = manifest["files"]
        if set(names) != set(expected) | {"reproducibility_manifest.json"}:
            raise ValueError("Archive contents do not match the manifest")
        for name, expected_hash in expected.items():
            digest = hashlib.sha256()
            with archive.open(name) as handle:
                for block in iter(lambda: handle.read(1024 * 1024), b""): digest.update(block)
            if digest.hexdigest() != expected_hash:
                raise ValueError(f"Checksum mismatch: {name}")
    return {"status": "VERIFIED", "file_count": len(expected), "scientific_accuracy_validated": False}
