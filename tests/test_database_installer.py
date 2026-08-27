from __future__ import annotations

import gzip
import hashlib
import io
import csv
import json
import tarfile
from pathlib import Path

import pytest

from phagemine import cli
from phagemine.database_installer import (
    DatabaseInstaller,
    RESOURCE_ORDER,
    installation_plan,
)
from phagemine.resources import EvidenceResourceManager, ResourceType


def _tar(path: Path, files: dict[str, bytes]) -> None:
    with tarfile.open(path, "w:gz") as archive:
        for name, content in files.items():
            info = tarfile.TarInfo(name)
            info.size = len(content)
            archive.addfile(info, io.BytesIO(content))


def test_cli_exposes_exact_install_commands(capsys):
    for resource in RESOURCE_ORDER:
        assert cli.main(["databases", "install", resource, "--dry-run"]) == 0
    assert cli.main(["databases", "install", "--all", "--dry-run"]) == 0
    text = capsys.readouterr().err
    assert "approximately" in text
    assert "3.4 GiB" in text


def test_no_argument_help_contains_post_install_guidance(capsys):
    assert cli.main([]) == 0
    output = capsys.readouterr().out
    assert "phagemine databases install --all" in output
    assert "phagemine doctor" in output


def test_installation_plan_lists_all_resources():
    plan = installation_plan(RESOURCE_ORDER)
    for resource in RESOURCE_ORDER:
        assert resource in plan
    assert "additional working and storage space" in plan


def test_pfam_install_prepares_registers_and_records_manifest(tmp_path, monkeypatch):
    version_file = tmp_path / "Pfam.version.gz"
    with gzip.open(version_file, "wt") as handle:
        handle.write("Pfam release       : 99.1\n")
    archive = tmp_path / "Pfam-A.hmm.gz"
    with gzip.open(archive, "wb") as handle:
        handle.write(b"HMMER3/f\nNAME  example\n//\n")
    md5 = hashlib.md5(archive.read_bytes()).hexdigest()
    checksums = tmp_path / "md5_checksums"
    checksums.write_text(f"{md5}  Pfam-A.hmm.gz\n")

    installer = DatabaseInstaller(
        tmp_path / "databases",
        registry_path=tmp_path / "registry.json",
        keep_downloads=True,
        reporter=lambda _message: None,
    )

    def fake_download(url, destination, **_kwargs):
        if url.endswith("Pfam.version.gz"):
            return version_file
        if url.endswith("md5_checksums"):
            return checksums
        return archive

    def fake_run(command):
        hmm = Path(command[-1])
        for suffix in (".h3f", ".h3i", ".h3m", ".h3p"):
            Path(str(hmm) + suffix).write_bytes(b"index")

    monkeypatch.setattr(installer, "_download", fake_download)
    monkeypatch.setattr(installer, "_require_tool", lambda name: name)
    monkeypatch.setattr(installer, "_run", fake_run)
    monkeypatch.setattr("phagemine.resources._tool_available", lambda _tool: True)

    result = installer.install("pfam")
    assert result.status == "READY"
    assert result.version == "99.1"
    assert Path(result.path).name == "Pfam-A.hmm"
    assert (Path(result.path).parent / "install_manifest.json").is_file()
    registered = installer.manager.validate("Pfam-A")
    assert registered["status"] == "READY"
    assert registered["provenance"]["installer"] == "phagemine databases install"


def test_vog_hmms_are_concatenated_in_pinned_archive_order(tmp_path):
    archive = tmp_path / "vog.tar.gz"
    _tar(archive, {"hmm/VOG0002.hmm": b"SECOND", "hmm/VOG0001.hmm": b"FIRST"})
    output = tmp_path / "VOGDB.hmm"
    DatabaseInstaller._concatenate_hmms(archive, output)
    assert output.read_bytes() == b"SECOND\nFIRST\n"


def test_phrogs_extraction_selects_only_required_safe_basenames(tmp_path):
    archive = tmp_path / "phrogs.tar.gz"
    _tar(archive, {
        "bundle/phrogs_profile_db": b"db",
        "bundle/phrogs_profile_db.dbtype": b"type",
        "bundle/phrogs_profile_db.index": b"index",
        "bundle/all_phrogs.h3m": b"hmm profiles",
        "bundle/phrog_annot_v4.tsv": b"phrog\tannot\n",
        "../../not_selected.txt": b"unsafe",
    })
    output = tmp_path / "output"
    output.mkdir()
    extracted = DatabaseInstaller._extract_phrogs(archive, output)
    assert extracted == [
        "all_phrogs.h3m", "phrog_annot_v4.tsv", "phrogs_profile_db",
        "phrogs_profile_db.dbtype", "phrogs_profile_db.index",
    ]
    assert not (tmp_path / "not_selected.txt").exists()


def test_uniprot_metalink_parser_reads_version_and_md5(tmp_path):
    path = tmp_path / "RELEASE.metalink"
    path.write_text("""<metalink xmlns="http://www.metalinker.org/" version="3.0">
      <version>2026_02</version><files>
      <file name="uniprot_sprot.fasta.gz"><verification><hash type="md5">abc</hash></verification></file>
      <file name="uniprot_sprot.dat.gz"><verification><hash type="md5">def</hash></verification></file>
      </files></metalink>""")
    version, checksums = DatabaseInstaller._parse_uniprot_metalink(path)
    assert version == "2026_02"
    assert checksums == {"uniprot_sprot.fasta.gz": "abc", "uniprot_sprot.dat.gz": "def"}


def test_cli_requires_resource_or_all():
    with pytest.raises(SystemExit):
        cli.main(["databases", "install", "--dry-run"])


def test_cli_attaches_existing_phrogs_hmm_without_download(tmp_path, monkeypatch):
    monkeypatch.setenv("PHAGEMINE_REGISTRY_PATH", str(tmp_path / "config" / "resources.json"))
    database = tmp_path / "phrogs_profile_db"
    database.write_bytes(b"db")
    Path(str(database) + ".dbtype").write_bytes(b"type")
    annotations = tmp_path / "phrog_annot_v4.tsv"
    annotations.write_text("phrog\tannot\n")
    hmm = tmp_path / "all_phrogs.h3m"
    hmm.write_bytes(b"profiles")
    manager = EvidenceResourceManager()
    manager.register(
        "PHROGs", ResourceType.PHROGS, database, version="v4",
        required_tools=["mmseqs"], preparation_status="prepared",
        provenance={"annotations_path": str(annotations)})
    monkeypatch.setattr("phagemine.resources._tool_available", lambda _tool: True)

    assert cli.main(["databases", "attach-phrogs-hmm", str(hmm)]) == 0

    registered = manager.validate("PHROGs")
    assert registered["status"] == "READY"
    assert registered["provenance"]["hmm_profiles_path"] == str(hmm.resolve())


def _gzip_text(path: Path, text: str) -> Path:
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        handle.write(text)
    return path


def test_pmfdb_install_converts_inphared_and_preserves_prediction_status(tmp_path, monkeypatch):
    proteins = _gzip_text(
        tmp_path / "proteins.faa.gz",
        ">A_00001 hypothetical protein\nMKT\n>A_00002 terminase large subunit\nMPEPTIDE\n",
    )
    mapping = _gzip_text(
        tmp_path / "mapping.csv.gz",
        "protein_id,contig_id,keywords\nA_00001,A,none\nA_00002,A,none\n",
    )
    metadata = _gzip_text(
        tmp_path / "metadata.tsv.gz",
        "Accession\tDescription\tGenome Length (KB)\tmolGC (%)\tGenus\tSub-family\tFamily\tHost\n"
        '<a href="https://example.test/A">A</a>\tExample phage\t40\t50\tTestvirus\tTestvirinae\tTestviridae\tPseudomonas\n',
    )
    installer = DatabaseInstaller(
        tmp_path / "databases", registry_path=tmp_path / "registry.json",
        keep_downloads=True, reporter=lambda _message: None,
    )

    def fake_download(url, destination, **_kwargs):
        if "proteins" in url:
            return proteins
        if "gene_to_genome" in url:
            return mapping
        return metadata

    def fake_run(command):
        if command[1] == "createdb":
            target = Path(command[3])
            target.write_bytes(b"db")
            Path(str(target) + ".dbtype").write_bytes(b"type")

    monkeypatch.setattr(installer, "_download", fake_download)
    monkeypatch.setattr(installer, "_require_tool", lambda name: name)
    monkeypatch.setattr(installer, "_run", fake_run)
    monkeypatch.setattr("phagemine.resources._tool_available", lambda _tool: True)
    result = installer.install("pmfdb")
    root = Path(result.path)
    assert result.status == "READY"
    records = list(csv.DictReader((root / "reference_metadata.tsv").open(), delimiter="\t"))
    assert records[0]["annotation_status"] == "PREDICTED_UNCHARACTERIZED"
    assert records[1]["annotation_status"] == "PREDICTED_FUNCTION"
    assert {row["characterized"] for row in records} == {"false"}
    manifest = json.loads((root / "reference_manifest.json").read_text())
    assert manifest["source_database"] == "INPHARED"
    assert manifest["protein_count"] == 2


def test_inphared_install_builds_genome_mash_resource(tmp_path, monkeypatch):
    genomes = _gzip_text(tmp_path / "genomes.fa.gz", ">A Example phage\nACGTACGT\n")
    metadata = _gzip_text(
        tmp_path / "metadata.tsv.gz",
        "Accession\tDescription\tGenome Length (KB)\tmolGC (%)\tGenus\tSub-family\tFamily\tHost\n"
        "A\tExample phage\t0.008\t50\tTestvirus\tTestvirinae\tTestviridae\tPseudomonas\n",
    )
    installer = DatabaseInstaller(
        tmp_path / "databases", registry_path=tmp_path / "registry.json",
        keep_downloads=True, reporter=lambda _message: None,
    )

    def fake_download(url, destination, **_kwargs):
        return genomes if "genomes.fa" in url else metadata

    def fake_run(command):
        if command[1] == "sketch":
            Path(command[command.index("-o") + 1] + ".msh").write_bytes(b"mash")

    monkeypatch.setattr(installer, "_download", fake_download)
    monkeypatch.setattr(installer, "_require_tool", lambda name: name)
    monkeypatch.setattr(installer, "_run", fake_run)
    monkeypatch.setattr("phagemine.resources._tool_available", lambda _tool: True)
    result = installer.install("inphared")
    root = Path(result.path).parent
    assert result.status == "READY"
    assert (root / "inphared.msh").is_file()
    manifest = json.loads((root / "genome_manifest.json").read_text())
    assert manifest["genome_count"] == 1
    install_manifest = json.loads((root / "install_manifest.json").read_text())
    assert "-i" in install_manifest["preparation_commands"][0]
