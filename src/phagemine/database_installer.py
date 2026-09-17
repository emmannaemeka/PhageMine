"""Download, prepare, register, and verify PhageMine evidence databases.

The installers deliberately use pinned or provider-described releases and keep
an installation manifest beside every prepared database.  Large biological
resources are not bundled with PhageMine.
"""
from __future__ import annotations

import gzip
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath
from typing import Callable, Iterable

from .inphared import INPHARED_RELEASE, build_genome_reference, build_pmfdb
from .resources import EvidenceResourceManager, ResourceType


PFAM_BASE = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release"
UNIPROT_BASE = "https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete"

# VOGDB is intentionally pinned. Updating this definition is an explicit,
# reviewable scientific change rather than a silent switch to a new release.
VOGDB_VERSION = "235"
VOGDB_BASE = f"https://fileshare.csb.univie.ac.at/vog/vog{VOGDB_VERSION}"
VOGDB_FILES = {
    "vog.hmm.tar.gz": "ef92a4c880126ea3b86e52fdf902ba06",
    "vog.annotations.tsv.gz": "1f6481d5b4fa82ec716caa142a638190",
}

# PHROGs v4 is distributed here as the MMseqs2-v18-compatible Pharokka 1.11.0
# database bundle. Only PHROGs files are extracted. The checksum is published
# by Pharokka's versioned database installer.
PHROGS_VERSION = "v4-pharokka-1.11.0"
PHROGS_ARCHIVE = "pharokka_v1.11.0_databases.tar.gz"
PHROGS_URL = f"https://zenodo.org/records/21755221/files/{PHROGS_ARCHIVE}"
PHROGS_MD5 = "143bb375ddb0b0653e5cb5671f4a7629"

INPHARED_BASE = "https://s3.climb.ac.uk/millardlab-inphared/2026"
INPHARED_FILES = {
    "7Apr2026_vConTACT2_proteins.faa.gz": {
        "md5": "e1551f1ffde752da8fc8d48d6ca7330c", "size": 329085310,
    },
    "7Apr2026_vConTACT2_gene_to_genome.csv.gz": {
        "md5": "32ab9c366d2c7e9580ae8fffc322fc70", "size": 5984549,
    },
    "7Apr2026_millardlab_website_table.txt.gz": {
        "md5": "7ca76b6d60f077a07d75e4c6ea635ea4", "size": 681798,
    },
    "7Apr2026_genomes.fa.gz": {
        "md5": "5fadc36048ced47f431fe510557cc797", "size": 714521461,
    },
}

RESOURCE_ORDER = ("pfam", "vogdb", "swissprot", "phrogs", "pmfdb", "inphared")
DOWNLOAD_ESTIMATES = {
    "pfam": 399 * 1024**2,
    "vogdb": 574 * 1024**2,
    "swissprot": 756 * 1024**2,
    "phrogs": 735 * 1024**2,
    "pmfdb": sum(INPHARED_FILES[name]["size"] for name in (
        "7Apr2026_vConTACT2_proteins.faa.gz",
        "7Apr2026_vConTACT2_gene_to_genome.csv.gz",
        "7Apr2026_millardlab_website_table.txt.gz",
    )),
    "inphared": sum(INPHARED_FILES[name]["size"] for name in (
        "7Apr2026_genomes.fa.gz",
        "7Apr2026_millardlab_website_table.txt.gz",
    )),
}


def default_database_root() -> Path:
    """Return a platform-appropriate per-user data directory."""
    system = platform.system()
    if system == "Darwin":
        return Path.home() / "Library" / "Application Support" / "PhageMine" / "databases"
    if system == "Windows":
        base = Path(os.environ.get("LOCALAPPDATA", Path.home() / "AppData" / "Local"))
        return base / "PhageMine" / "databases"
    base = Path(os.environ.get("XDG_DATA_HOME", Path.home() / ".local" / "share"))
    return base / "phagemine" / "databases"


def human_size(value: int) -> str:
    amount = float(value)
    for unit in ("B", "KiB", "MiB", "GiB", "TiB"):
        if amount < 1024 or unit == "TiB":
            return f"{amount:.1f} {unit}"
        amount /= 1024
    return f"{amount:.1f} TiB"


def installation_plan(resources: Iterable[str]) -> str:
    names = tuple(resources)
    total = sum(DOWNLOAD_ESTIMATES[name] for name in names)
    lines = ["PhageMine database installation plan:"]
    lines.extend(f"  {name:<11} approximately {human_size(DOWNLOAD_ESTIMATES[name])} compressed" for name in names)
    lines.append(f"  {'total':<11} approximately {human_size(total)} compressed")
    lines.append("Prepared databases require additional working and storage space.")
    return "\n".join(lines)


@dataclass(frozen=True)
class InstallResult:
    resource: str
    version: str
    path: str
    status: str
    installed: bool


class DatabaseInstallError(RuntimeError):
    """Raised when a database cannot be installed safely."""


class DatabaseInstaller:
    def __init__(
        self,
        directory: str | Path | None = None,
        *,
        registry_path: str | Path | None = None,
        force: bool = False,
        keep_downloads: bool = False,
        threads: int = 1,
        reporter: Callable[[str], None] | None = None,
    ):
        self.directory = Path(directory).expanduser().resolve() if directory else default_database_root()
        self.manager = EvidenceResourceManager(registry_path)
        self.force = force
        self.keep_downloads = keep_downloads
        self.threads = max(1, int(threads))
        self.reporter = reporter or (lambda message: print(message, file=sys.stderr, flush=True))
        self.download_dir = self.directory / ".downloads"

    def install_many(self, resources: Iterable[str]) -> list[InstallResult]:
        return [self.install(name) for name in resources]

    def install(self, resource: str) -> InstallResult:
        name = resource.lower()
        installers = {
            "pfam": self._install_pfam,
            "vogdb": self._install_vogdb,
            "swissprot": self._install_swissprot,
            "phrogs": self._install_phrogs,
            "pmfdb": self._install_pmfdb,
            "inphared": self._install_inphared,
        }
        if name not in installers:
            raise DatabaseInstallError(f"unsupported database: {resource}")
        self.directory.mkdir(parents=True, exist_ok=True)
        self.download_dir.mkdir(parents=True, exist_ok=True)
        return installers[name]()

    def _already_ready(self, registry_name: str, resource: str | None = None) -> InstallResult | None:
        current = self.manager.validate(registry_name)
        if current and current.get("status") == "READY" and not self.force:
            self.reporter(f"{registry_name} is already READY at {current['path']}")
            return InstallResult(resource or registry_name, current.get("version") or "unknown", current["path"], "READY", False)
        return None

    def _transaction(
        self,
        *,
        resource: str,
        registry_name: str,
        resource_type: ResourceType,
        version: str,
        required_tools: list[str],
        prepare: Callable[[Path], dict],
    ) -> InstallResult:
        ready = self._already_ready(registry_name, resource)
        if ready:
            return ready
        resource_root = self.directory / resource
        target = resource_root / version
        resource_root.mkdir(parents=True, exist_ok=True)
        if target.exists() and not self.force:
            raise DatabaseInstallError(
                f"installation directory already exists but is not registered READY: {target}; "
                "inspect it or rerun with --force"
            )
        stage = Path(tempfile.mkdtemp(prefix=".install-", dir=resource_root))
        try:
            details = prepare(stage)
            manifest = {
                "format_version": "1",
                "resource": resource,
                "resource_type": resource_type.value,
                "version": version,
                "installed_at": datetime.now(timezone.utc).isoformat(),
                **details,
            }
            (stage / "install_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
            primary_relative = Path(details["primary_path"])
            if not (stage / primary_relative).exists():
                raise DatabaseInstallError(f"installer did not create expected primary database: {primary_relative}")
            if target.exists():
                self._remove_target(target, resource_root)
            stage.replace(target)
        except Exception:
            shutil.rmtree(stage, ignore_errors=True)
            raise

        provenance = dict(details.get("provenance") or {})
        for key, value in list(provenance.items()):
            if key.endswith("_path") and value and not Path(str(value)).is_absolute():
                provenance[key] = str(target / str(value))
        provenance["install_manifest"] = str(target / "install_manifest.json")
        provenance["installer"] = "phagemine databases install"
        primary = target / details["primary_path"]
        registered = self.manager.register(
            registry_name,
            resource_type,
            primary,
            version=version,
            required_tools=required_tools,
            preparation_status="prepared",
            notes="Installed and prepared by PhageMine",
            provenance=provenance,
        )
        checked = self.manager.validate(registry_name)
        if not checked or checked.get("status") != "READY":
            errors = ", ".join((checked or {}).get("validation_errors", [])) or "unknown validation error"
            self.manager.unregister(registry_name)
            raise DatabaseInstallError(f"{registry_name} was installed but failed validation: {errors}")
        self.reporter(f"{registry_name} {version} installed and registered: {registered.path}")
        return InstallResult(resource, version, registered.path, "READY", True)

    @staticmethod
    def _remove_target(target: Path, resource_root: Path) -> None:
        resolved_target = target.resolve()
        resolved_root = resource_root.resolve()
        if resolved_target.parent != resolved_root:
            raise DatabaseInstallError(f"refusing to replace unexpected path: {target}")
        if target.is_dir():
            shutil.rmtree(target)
        else:
            target.unlink()

    def _require_tool(self, name: str) -> str:
        found = shutil.which(name)
        if not found:
            raise DatabaseInstallError(
                f"required executable '{name}' is unavailable. Install the PhageMine Conda dependencies and retry."
            )
        return found

    def _run(self, command: list[str]) -> None:
        self.reporter("Preparing database: " + " ".join(command))
        try:
            subprocess.run(command, check=True, text=True, capture_output=True)
        except subprocess.CalledProcessError as exc:
            detail = (exc.stderr or exc.stdout or "").strip()
            raise DatabaseInstallError(f"database preparation failed ({' '.join(command)}): {detail}") from exc

    def _download(
        self,
        url: str,
        destination: Path,
        *,
        expected: str | None = None,
        algorithm: str = "md5",
        refresh: bool = False,
    ) -> Path:
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() and not refresh:
            if not expected or self._digest(destination, algorithm) == expected.lower():
                self.reporter(f"Using cached download: {destination.name}")
                return destination
            destination.unlink()
        part = destination.with_name(destination.name + ".part")
        offset = part.stat().st_size if part.exists() else 0
        headers = {"User-Agent": "PhageMine/1.0 database-installer"}
        if offset:
            headers["Range"] = f"bytes={offset}-"
        self.reporter(f"Downloading {url} (connection/read timeout: 60 seconds)")
        try:
            request = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(request, timeout=60) as response:
                status = getattr(response, "status", None) or response.getcode()
                append = bool(offset and status == 206)
                if not append:
                    offset = 0
                total_header = response.headers.get("Content-Length")
                total = (int(total_header) + offset) if total_header else None
                mode = "ab" if append else "wb"
                downloaded = offset
                started = last_report = time.monotonic()

                def report_progress(label: str) -> None:
                    elapsed = max(time.monotonic() - started, 0.001)
                    speed = (downloaded - offset) / elapsed
                    detail = human_size(downloaded)
                    if total:
                        detail += f" / {human_size(total)} ({100 * downloaded / total:.1f}%)"
                    else:
                        detail += " (total size unknown)"
                    detail += f" | {human_size(speed)}/s"
                    if total and speed > 0:
                        detail += f" | ETA {max(0, (total - downloaded) / speed):.0f}s"
                    self.reporter(f"  {label}: {detail}")

                report_progress("Resuming" if append else "Starting")
                # read1 returns available data without waiting to fill a MiB buffer.
                read_chunk = getattr(response, "read1", response.read)
                with part.open(mode) as handle:
                    while True:
                        chunk = read_chunk(1024 * 1024)
                        if not chunk:
                            break
                        handle.write(chunk)
                        downloaded += len(chunk)
                        now = time.monotonic()
                        if now - last_report >= 2:
                            report_progress("Received")
                            last_report = now
                if total is not None and downloaded != total:
                    raise OSError(f"incomplete download: received {downloaded} of {total} bytes")
                report_progress("Download complete")
        except (OSError, urllib.error.URLError) as exc:
            raise DatabaseInstallError(f"download failed; partial data retained for resume: {url}: {exc}") from exc
        part.replace(destination)
        if expected:
            self.reporter(f"Verifying {algorithm} checksum: {destination.name}")
            observed = self._digest(destination, algorithm)
            if observed != expected.lower():
                destination.unlink(missing_ok=True)
                raise DatabaseInstallError(
                    f"checksum mismatch for {destination.name}: expected {expected.lower()}, observed {observed}"
                )
        return destination

    @staticmethod
    def _digest(path: Path, algorithm: str) -> str:
        digest = hashlib.new(algorithm)
        with path.open("rb") as handle:
            for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(chunk)
        return digest.hexdigest()

    @staticmethod
    def _gunzip(source: Path, destination: Path) -> None:
        with gzip.open(source, "rb") as input_handle, destination.open("wb") as output_handle:
            shutil.copyfileobj(input_handle, output_handle, length=1024 * 1024)

    def _cleanup(self, paths: Iterable[Path]) -> None:
        if self.keep_downloads:
            return
        for path in paths:
            path.unlink(missing_ok=True)

    def _install_pfam(self) -> InstallResult:
        ready = self._already_ready("Pfam-A", "pfam")
        if ready:
            return ready
        version_file = self._download(
            f"{PFAM_BASE}/Pfam.version.gz", self.download_dir / "Pfam.version.gz", refresh=True
        )
        checksum_file = self._download(
            f"{PFAM_BASE}/md5_checksums", self.download_dir / "Pfam.md5_checksums", refresh=True
        )
        with gzip.open(version_file, "rt", encoding="utf-8") as handle:
            version_text = handle.read()
        version = next(
            (line.split(":", 1)[1].strip() for line in version_text.splitlines() if line.startswith("Pfam release")),
            None,
        )
        if not version:
            raise DatabaseInstallError("could not determine Pfam release from Pfam.version.gz")
        checksums = self._parse_checksum_file(checksum_file.read_text())
        expected = checksums.get("Pfam-A.hmm.gz")
        if not expected:
            raise DatabaseInstallError("official Pfam checksum list does not contain Pfam-A.hmm.gz")
        archive = self._download(
            f"{PFAM_BASE}/Pfam-A.hmm.gz",
            self.download_dir / f"Pfam-{version}-A.hmm.gz",
            expected=expected,
        )
        hmmpress = self._require_tool("hmmpress")

        def prepare(stage: Path) -> dict:
            hmm = stage / "Pfam-A.hmm"
            self._gunzip(archive, hmm)
            self._run([hmmpress, "-f", str(hmm)])
            self._require_hmm_indexes(hmm)
            (stage / "Pfam.version.txt").write_text(version_text)
            shutil.copy2(checksum_file, stage / "md5_checksums")
            return {
                "primary_path": hmm.name,
                "source_artifacts": [{"url": f"{PFAM_BASE}/Pfam-A.hmm.gz", "md5": expected}],
                "preparation_commands": [["hmmpress", "-f", hmm.name]],
                "provenance": {"provider": "Pfam/InterPro", "release": version},
            }

        result = self._transaction(
            resource="pfam", registry_name="Pfam-A", resource_type=ResourceType.PFAM,
            version=version, required_tools=["hmmscan"], prepare=prepare,
        )
        self._cleanup([archive])
        return result

    def _install_vogdb(self) -> InstallResult:
        ready = self._already_ready("VOGDB", "vogdb")
        if ready:
            return ready
        archive = self._download(
            f"{VOGDB_BASE}/vog.hmm.tar.gz",
            self.download_dir / f"vog{VOGDB_VERSION}.hmm.tar.gz",
            expected=VOGDB_FILES["vog.hmm.tar.gz"],
        )
        annotations = self._download(
            f"{VOGDB_BASE}/vog.annotations.tsv.gz",
            self.download_dir / f"vog{VOGDB_VERSION}.annotations.tsv.gz",
            expected=VOGDB_FILES["vog.annotations.tsv.gz"],
        )
        hmmpress = self._require_tool("hmmpress")

        def prepare(stage: Path) -> dict:
            hmm = stage / "VOGDB.hmm"
            self._concatenate_hmms(archive, hmm)
            self._run([hmmpress, "-f", str(hmm)])
            self._require_hmm_indexes(hmm)
            shutil.copy2(annotations, stage / "vog.annotations.tsv.gz")
            return {
                "primary_path": hmm.name,
                "source_artifacts": [
                    {"url": f"{VOGDB_BASE}/vog.hmm.tar.gz", "md5": VOGDB_FILES["vog.hmm.tar.gz"]},
                    {"url": f"{VOGDB_BASE}/vog.annotations.tsv.gz", "md5": VOGDB_FILES["vog.annotations.tsv.gz"]},
                ],
                "preparation_commands": [["concatenate", "*.hmm", hmm.name], ["hmmpress", "-f", hmm.name]],
                "provenance": {
                    "provider": "VOGDB", "release": VOGDB_VERSION,
                    "annotations_path": "vog.annotations.tsv.gz",
                },
            }

        result = self._transaction(
            resource="vogdb", registry_name="VOGDB", resource_type=ResourceType.VOGDB,
            version=VOGDB_VERSION, required_tools=["hmmscan"], prepare=prepare,
        )
        self._cleanup([archive, annotations])
        return result

    def _install_swissprot(self) -> InstallResult:
        ready = self._already_ready("Swiss-Prot", "swissprot")
        if ready:
            return ready
        metalink = self._download(
            f"{UNIPROT_BASE}/RELEASE.metalink",
            self.download_dir / "UniProt-RELEASE.metalink",
            refresh=True,
        )
        version, checksums = self._parse_uniprot_metalink(metalink)
        required = ("uniprot_sprot.fasta.gz", "uniprot_sprot.dat.gz")
        if any(name not in checksums for name in required):
            raise DatabaseInstallError("UniProt RELEASE.metalink is missing Swiss-Prot checksums")
        fasta = self._download(
            f"{UNIPROT_BASE}/uniprot_sprot.fasta.gz",
            self.download_dir / f"uniprot_sprot-{version}.fasta.gz",
            expected=checksums["uniprot_sprot.fasta.gz"],
        )
        metadata = self._download(
            f"{UNIPROT_BASE}/uniprot_sprot.dat.gz",
            self.download_dir / f"uniprot_sprot-{version}.dat.gz",
            expected=checksums["uniprot_sprot.dat.gz"],
        )
        diamond = self._require_tool("diamond")

        def prepare(stage: Path) -> dict:
            prefix = stage / "uniprot_sprot"
            self._run([diamond, "makedb", "--in", str(fasta), "--db", str(prefix)])
            dmnd = stage / "uniprot_sprot.dmnd"
            if not dmnd.is_file():
                raise DatabaseInstallError("DIAMOND did not create uniprot_sprot.dmnd")
            shutil.copy2(metadata, stage / "uniprot_sprot.dat.gz")
            shutil.copy2(metalink, stage / "RELEASE.metalink")
            return {
                "primary_path": dmnd.name,
                "source_artifacts": [
                    {"url": f"{UNIPROT_BASE}/{name}", "md5": checksums[name]} for name in required
                ],
                "preparation_commands": [["diamond", "makedb", "--in", fasta.name, "--db", "uniprot_sprot"]],
                "provenance": {
                    "provider": "UniProtKB/Swiss-Prot", "release": version,
                    "metadata_path": "uniprot_sprot.dat.gz",
                },
            }

        result = self._transaction(
            resource="swissprot", registry_name="Swiss-Prot", resource_type=ResourceType.SWISSPROT,
            version=version, required_tools=["diamond"], prepare=prepare,
        )
        self._cleanup([fasta, metadata])
        return result

    def _install_phrogs(self) -> InstallResult:
        ready = self._already_ready("PHROGs", "phrogs")
        if ready:
            return ready
        self._require_tool("mmseqs")
        archive = self._download(
            PHROGS_URL, self.download_dir / PHROGS_ARCHIVE, expected=PHROGS_MD5
        )

        def prepare(stage: Path) -> dict:
            extracted = self._extract_phrogs(archive, stage)
            primary = stage / "phrogs_profile_db"
            hmm_profiles = stage / "all_phrogs.h3m"
            annotations = stage / "phrog_annot_v4.tsv"
            if not primary.is_file() or not Path(str(primary) + ".dbtype").is_file():
                raise DatabaseInstallError("PHROGs bundle does not contain a prepared MMseqs2 profile database")
            if not annotations.is_file():
                raise DatabaseInstallError("PHROGs bundle does not contain phrog_annot_v4.tsv")
            if not hmm_profiles.is_file():
                raise DatabaseInstallError("PHROGs bundle does not contain all_phrogs.h3m")
            return {
                "primary_path": primary.name,
                "source_artifacts": [{"url": PHROGS_URL, "md5": PHROGS_MD5}],
                "extracted_files": extracted,
                "preparation_commands": [],
                "provenance": {
                    "provider": "PHROGs", "release": "v4",
                    "distribution": "Pharokka database 1.11.0",
                    "annotations_path": annotations.name,
                    "hmm_profiles_path": hmm_profiles.name,
                },
            }

        result = self._transaction(
            resource="phrogs", registry_name="PHROGs", resource_type=ResourceType.PHROGS,
            version=PHROGS_VERSION, required_tools=["mmseqs"], prepare=prepare,
        )
        self._cleanup([archive])
        return result

    @staticmethod
    def _inphared_artifact(name: str) -> dict:
        definition = INPHARED_FILES[name]
        return {"url": f"{INPHARED_BASE}/{name}", "md5": definition["md5"], "bytes": definition["size"]}

    def _inphared_download(self, name: str) -> Path:
        return self._download(
            f"{INPHARED_BASE}/{name}",
            self.download_dir / name,
            expected=INPHARED_FILES[name]["md5"],
        )

    def _install_pmfdb(self) -> InstallResult:
        ready = self._already_ready("PMFDB-INPHARED", "pmfdb")
        if ready:
            return ready
        protein_name = "7Apr2026_vConTACT2_proteins.faa.gz"
        mapping_name = "7Apr2026_vConTACT2_gene_to_genome.csv.gz"
        metadata_name = "7Apr2026_millardlab_website_table.txt.gz"
        proteins = self._inphared_download(protein_name)
        mapping = self._inphared_download(mapping_name)
        metadata = self._inphared_download(metadata_name)
        mmseqs = self._require_tool("mmseqs")
        artifacts = [self._inphared_artifact(name) for name in (protein_name, mapping_name, metadata_name)]

        def prepare(stage: Path) -> dict:
            manifest = build_pmfdb(
                proteins, mapping, metadata, stage,
                mmseqs=mmseqs, run=self._run, threads=self.threads,
                source_artifacts=artifacts, release=INPHARED_RELEASE,
            )
            return {
                "primary_path": ".",
                "source_artifacts": artifacts,
                "preparation_commands": [
                    ["mmseqs", "createdb", "reference_phage_proteins.faa", "mmseqs/target_db"],
                    ["mmseqs", "createindex", "mmseqs/target_db", "mmseqs/tmp", "--threads", str(self.threads)],
                ],
                "provenance": {
                    "provider": "INPHARED",
                    "release": INPHARED_RELEASE,
                    "pmfdb_version": manifest["pmfdb_version"],
                    "metadata_path": "reference_metadata.tsv",
                    "qc_path": "reference_qc.tsv",
                    "reference_manifest_path": "reference_manifest.json",
                    "mmseqs_target_path": "mmseqs/target_db",
                },
            }

        result = self._transaction(
            resource="pmfdb", registry_name="PMFDB-INPHARED", resource_type=ResourceType.PMFDB,
            version=INPHARED_RELEASE, required_tools=["mmseqs"], prepare=prepare,
        )
        self._cleanup([proteins, mapping, metadata])
        return result

    def _install_inphared(self) -> InstallResult:
        ready = self._already_ready("INPHARED-Genomes", "inphared")
        if ready:
            return ready
        genome_name = "7Apr2026_genomes.fa.gz"
        metadata_name = "7Apr2026_millardlab_website_table.txt.gz"
        genomes = self._inphared_download(genome_name)
        metadata = self._inphared_download(metadata_name)
        mash = self._require_tool("mash")
        artifacts = [self._inphared_artifact(name) for name in (genome_name, metadata_name)]

        def prepare(stage: Path) -> dict:
            manifest = build_genome_reference(
                genomes, metadata, stage, mash=mash, run=self._run,
                source_artifacts=artifacts, release=INPHARED_RELEASE,
            )
            return {
                "primary_path": "reference_phage_genomes.fna",
                "source_artifacts": artifacts,
                "preparation_commands": [["mash", "sketch", "-i", "-o", "inphared", "reference_phage_genomes.fna"]],
                "provenance": {
                    "provider": "INPHARED",
                    "release": INPHARED_RELEASE,
                    "database_version": manifest["database_version"],
                    "metadata_path": "genome_metadata.tsv",
                    "qc_path": "genome_qc.tsv",
                    "reference_manifest_path": "genome_manifest.json",
                    "mash_index_path": "inphared.msh",
                },
            }

        result = self._transaction(
            resource="inphared", registry_name="INPHARED-Genomes",
            resource_type=ResourceType.INPHARED_GENOMES,
            version=INPHARED_RELEASE, required_tools=["mash"], prepare=prepare,
        )
        self._cleanup([genomes, metadata])
        return result

    @staticmethod
    def _parse_checksum_file(text: str) -> dict[str, str]:
        result: dict[str, str] = {}
        for line in text.splitlines():
            fields = line.split()
            if len(fields) >= 2:
                result[fields[-1].lstrip("*")] = fields[0].lower()
        return result

    @staticmethod
    def _parse_uniprot_metalink(path: Path) -> tuple[str, dict[str, str]]:
        root = ET.parse(path).getroot()
        version_node = root.find("{*}version")
        if version_node is None or not version_node.text:
            raise DatabaseInstallError("could not determine UniProt release")
        checksums: dict[str, str] = {}
        for file_node in root.findall(".//{*}file"):
            name = file_node.attrib.get("name")
            for hash_node in file_node.findall(".//{*}hash"):
                if name and hash_node.attrib.get("type", "").lower() == "md5" and hash_node.text:
                    checksums[name] = hash_node.text.strip().lower()
        return version_node.text.strip(), checksums

    @staticmethod
    def _require_hmm_indexes(hmm: Path) -> None:
        missing = [str(hmm) + suffix for suffix in (".h3f", ".h3i", ".h3m", ".h3p")
                   if not Path(str(hmm) + suffix).is_file()]
        if missing:
            raise DatabaseInstallError("hmmpress did not create all required indexes: " + ", ".join(missing))

    @staticmethod
    def _concatenate_hmms(archive: Path, destination: Path) -> None:
        with tarfile.open(archive, "r:gz") as tar, destination.open("wb") as output:
            count = 0
            # Preserve the order of the checksum-pinned archive. Streaming in
            # archive order avoids pathological repeated seeks through a large
            # gzip file while remaining reproducible for this exact snapshot.
            for member in tar:
                if not member.isfile() or not member.name.endswith(".hmm"):
                    continue
                source = tar.extractfile(member)
                if source is None:
                    raise DatabaseInstallError(f"could not read {member.name} from VOGDB archive")
                shutil.copyfileobj(source, output, length=1024 * 1024)
                output.write(b"\n")
                count += 1
            if not count:
                raise DatabaseInstallError("VOGDB archive contains no HMM files")

    @staticmethod
    def _extract_phrogs(archive: Path, destination: Path) -> list[str]:
        extracted: list[str] = []
        with tarfile.open(archive, "r:gz") as tar:
            for member in tar:
                if not member.isfile():
                    continue
                basename = PurePosixPath(member.name).name
                wanted = (
                    basename == "phrog_annot_v4.tsv"
                    or basename == "all_phrogs.h3m"
                    or basename.startswith("phrogs_profile_db")
                    or basename == "VERSION_1_11_0"
                )
                if not wanted:
                    continue
                target = destination / basename
                if target.exists():
                    raise DatabaseInstallError(f"duplicate PHROGs bundle member: {basename}")
                source = tar.extractfile(member)
                if source is None:
                    raise DatabaseInstallError(f"could not read {member.name} from PHROGs archive")
                with target.open("wb") as output:
                    shutil.copyfileobj(source, output, length=1024 * 1024)
                extracted.append(basename)
        return sorted(extracted)


def results_json(results: Iterable[InstallResult]) -> str:
    return json.dumps([asdict(result) for result in results], indent=2, sort_keys=True)
