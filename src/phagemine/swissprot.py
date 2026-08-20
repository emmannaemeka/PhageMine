"""Optional local Swiss-Prot curated sequence-homology evidence adapter."""
from __future__ import annotations

import gzip, json, sqlite3, hashlib, os
import re
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from .evidence import EvidenceAdapter, EvidenceAdapterResult
from .models import Evidence, EvidenceLevel, Protein


class SwissProtEvidenceAdapter(EvidenceAdapter):
    name = "SwissProtEvidenceAdapter"

    def __init__(self, database_path: str | Path | None = None, metadata_path: str | Path | None = None,
                 diamond: str | None = None, database_version: str | None = None,
                 strong_identity: float = 0.70, moderate_identity: float = 0.35,
                 strong_coverage: float = 0.70, moderate_coverage: float = 0.40,
                 evalue_threshold: float = 1e-5, min_alignment_length: int = 30, threads: int = 1):
        self.database_path = Path(database_path).expanduser() if database_path else None
        self.metadata_path = Path(metadata_path).expanduser() if metadata_path else None
        self.diamond = diamond or shutil.which("diamond")
        self.database_version = database_version or "unknown"
        self.strong_identity, self.moderate_identity = strong_identity, moderate_identity
        self.strong_coverage, self.moderate_coverage = strong_coverage, moderate_coverage
        self.evalue_threshold = evalue_threshold
        self.min_alignment_length = min_alignment_length
        self.threads = max(1, int(threads))

    def available(self) -> bool:
        executable = bool(self.diamond and (Path(self.diamond).exists() or shutil.which(self.diamond)))
        return bool(executable and self.database_path and self.database_path.is_file())

    def version(self) -> str:
        if not self.diamond or not (Path(self.diamond).exists() or shutil.which(self.diamond)):
            return "unavailable"
        result = subprocess.run([self.diamond, "version"], capture_output=True, text=True, check=False)
        match = re.search(r"([0-9]+\.[0-9]+(?:\.[0-9]+)?)", (result.stdout or result.stderr))
        return match.group(1) if match else "unknown"

    def provenance(self) -> dict[str, Any]:
        return {"adapter": self.name, "adapter_version": "0.1.0", "diamond": self.diamond,
                "diamond_version": self.version(), "swissprot_path": str(self.database_path) if self.database_path else None,
                "metadata_path": str(self.metadata_path) if self.metadata_path else None,
                "swissprot_version": self.database_version, "status": "REAL" if self.available() else "UNAVAILABLE",
                "thresholds": {"strong_identity": self.strong_identity, "moderate_identity": self.moderate_identity,
                               "strong_coverage": self.strong_coverage, "moderate_coverage": self.moderate_coverage,
                               "evalue": self.evalue_threshold, "min_alignment_length": self.min_alignment_length}}

    @staticmethod
    def parse_subject_id(subject_id: str) -> tuple[str | None, str | None]:
        parts = subject_id.split("|")
        if len(parts) >= 3 and parts[0] == "sp":
            return parts[1], parts[2]
        return subject_id, None

    def analyze(self, proteins: list[Protein]) -> EvidenceAdapterResult:
        provenance = self.provenance()
        if not self.available():
            return EvidenceAdapterResult(self.name, "UNAVAILABLE", provenance=provenance,
                message="Swiss-Prot/DIAMOND unavailable; no curated homology evidence was fabricated.")
        with tempfile.TemporaryDirectory(prefix="phagemine-swissprot-") as temp:
            fasta = Path(temp) / "proteins.faa"
            output = Path(temp) / "diamond.tsv"
            fasta.write_text("".join(f">{p.protein_id}\n{p.sequence}\n" for p in proteins))
            fields = "qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore"
            command = [str(self.diamond), "blastp", "--db", str(self.database_path), "--query", str(fasta),
                       "--out", str(output), "--outfmt", "6", "--threads", str(self.threads), *fields.split()]
            result = subprocess.run(command, capture_output=True, text=True, check=False)
            provenance["search_command"] = command
            provenance["search_parameters"] = {"outfmt": fields}
            if result.returncode:
                return EvidenceAdapterResult(self.name, "UNAVAILABLE", provenance=provenance,
                    message=f"DIAMOND failed with exit code {result.returncode}: {(result.stderr or result.stdout).strip()}")
            hits = self.parse_tabular(output.read_text(), proteins, provenance)
        self._mark_conflicts(hits)
        return EvidenceAdapterResult(self.name, "REAL", evidence=hits, provenance=provenance)

    def _classify(self, identity: float, qcov: float, scov: float, evalue: float, length: int) -> str:
        if evalue > self.evalue_threshold or length < self.min_alignment_length:
            return "REJECTED"
        if identity >= self.strong_identity and qcov >= self.strong_coverage and scov >= self.strong_coverage:
            return "STRONG"
        if identity >= self.moderate_identity and qcov >= self.moderate_coverage and scov >= self.moderate_coverage:
            return "MODERATE"
        return "REJECTED"

    def parse_tabular(self, text: str, proteins: list[Protein], provenance: dict[str, Any] | None = None) -> list[Evidence]:
        lengths = {p.protein_id: p.length for p in proteins}
        records: list[Evidence] = []
        for line in text.splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 12:
                continue
            try:
                protein_id, subject = fields[0], fields[1]
                identity, alignment_length = float(fields[2]), int(fields[3])
                query_length, subject_length = int(fields[4]), int(fields[5])
                qstart, qend, sstart, send = (int(fields[i]) for i in range(6, 10))
                evalue, bitscore = float(fields[10]), float(fields[11])
            except (ValueError, IndexError):
                continue
            if protein_id not in lengths:
                continue
            accession, entry_name = self.parse_subject_id(subject)
            # DIAMOND's alignment length includes gapped columns; coverage is
            # based on the inclusive aligned coordinate spans instead.
            query_span = abs(qend - qstart) + 1
            subject_span = abs(send - sstart) + 1
            qcov = query_span / query_length if query_length else 0.0
            scov = subject_span / subject_length if subject_length else 0.0
            qcov = min(1.0, max(0.0, qcov))
            scov = min(1.0, max(0.0, scov))
            strength = self._classify(identity / 100.0, qcov, scov, evalue, alignment_length)
            metadata = self._metadata_for(accession)
            supports = strength in {"STRONG", "MODERATE"}
            metrics = {"accession": accession, "entry_name": entry_name, "reviewed": metadata.get("reviewed"),
                       "protein_name": metadata.get("protein_name"), "organism": metadata.get("organism"),
                       "taxonomy_id": metadata.get("taxonomy_id"), "gene_name": metadata.get("gene_name"),
                       "protein_existence": metadata.get("protein_existence"), "percent_identity": identity,
                       "alignment_length": alignment_length, "query_length": query_length, "subject_length": subject_length,
                       "query_coordinates": {"start": qstart, "end": qend}, "subject_coordinates": {"start": sstart, "end": send},
                       "query_coverage": qcov, "subject_coverage": scov, "evalue": evalue, "e_value": evalue, "bit_score": bitscore,
                       "raw_diamond_row": line, "conflict": False}
            records.append(Evidence("sequence_similarity", f"Swiss-Prot curated homology evidence for {accession}; evidence only, not an automatic product assignment.",
                EvidenceLevel.CURATED if strength == "STRONG" else EvidenceLevel.COMPUTATIONAL if strength == "MODERATE" else EvidenceLevel.WEAK,
                "Swiss-Prot", self.database_version, status="REAL", supports=supports, metrics=metrics, identifier=accession,
                threshold={"evalue": self.evalue_threshold, "strong_identity": self.strong_identity, "moderate_identity": self.moderate_identity},
                coordinates={"start": qstart, "end": qend}, provenance={**(provenance or {}), "protein_id": protein_id, "accession": accession},
                evidence_strength=strength, family_name=entry_name, description=metadata.get("protein_name")))
        return records

    def _metadata_for(self, accession: str | None) -> dict[str, Any]:
        if not accession or not self.metadata_path or not self.metadata_path.exists():
            return {}
        index=self._metadata_index()
        if index:
            with sqlite3.connect(index) as db:
                row=db.execute('SELECT payload FROM metadata WHERE accession=?',(accession,)).fetchone()
            return json.loads(row[0]) if row else {}
        opener = gzip.open if self.metadata_path.suffix == ".gz" else open
        current: list[str] = []
        result: dict[str, Any] = {}
        with opener(self.metadata_path, "rt", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("//"):
                    if current and any(x == f"AC   {accession};" for x in current):
                        result = self._parse_dat_record(current)
                        break
                    current = []
                current.append(line.rstrip("\n"))
        return result

    def _metadata_index(self) -> Path | None:
        """Build/reuse a checksum-validated local metadata index atomically."""
        if not self.metadata_path or not self.metadata_path.exists(): return None
        source=str(self.metadata_path); stat=self.metadata_path.stat(); sha=hashlib.sha256(self.metadata_path.read_bytes()).hexdigest()
        index=self.metadata_path.with_suffix(self.metadata_path.suffix+'.sqlite')
        try:
            db = sqlite3.connect(index)
            try:
                meta=db.execute("SELECT value FROM meta WHERE key='source_sha256'").fetchone()
                if meta and meta[0]==sha: return index
            finally:
                db.close()
        except sqlite3.Error: pass
        tmp=index.with_name(index.name+'.tmp')
        if tmp.exists(): tmp.unlink()
        db = sqlite3.connect(tmp)
        try:
            db.execute('CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT)'); db.execute('CREATE TABLE metadata (accession TEXT PRIMARY KEY, payload TEXT)')
            opener=gzip.open if self.metadata_path.suffix=='.gz' else open; current=[]; rows=[]
            with opener(self.metadata_path,'rt',encoding='utf-8',errors='replace') as handle:
                for line in handle:
                    if line.startswith('//'):
                        if current:
                            rec=self._parse_dat_record(current)
                            ac=next((x[5:].split(';')[0] for x in current if x.startswith('AC   ')),None)
                            if ac: rows.append((ac,json.dumps(rec,sort_keys=True)))
                        current=[]
                    else: current.append(line.rstrip('\n'))
            db.executemany('INSERT OR REPLACE INTO metadata VALUES (?,?)',rows)
            db.executemany('INSERT INTO meta VALUES (?,?)',[('source_sha256',sha),('source_path',source),('source_size',str(stat.st_size)),('index_version','1')]); db.commit()
        finally:
            db.close()
        os.replace(tmp,index); return index

    @staticmethod
    def _parse_dat_record(lines: list[str]) -> dict[str, Any]:
        result: dict[str, Any] = {"reviewed": False}
        for line in lines:
            if line.startswith("ID   "):
                result["reviewed"] = "Reviewed" in line
            elif line.startswith("DE   RecName: Full="):
                result["protein_name"] = line.split("Full=", 1)[1].rstrip(";")
            elif line.startswith("OS   ") and "organism" not in result:
                result["organism"] = line[5:].strip()
            elif line.startswith("OX   NCBI_TaxID="):
                result["taxonomy_id"] = line.split("=", 1)[1].split(";", 1)[0]
            elif line.startswith("GN   "):
                match = re.search(r"(?:Name|ORFNames)=([^; ,]+)", line)
                if match:
                    result["gene_name"] = match.group(1)
            elif line.startswith("PE   "):
                result["protein_existence"] = line[5:].strip().rstrip(";")
        return result

    @staticmethod
    def _mark_conflicts(records: list[Evidence]) -> None:
        by_protein: dict[str, set[str]] = {}
        for record in records:
            if record.supports and record.description:
                by_protein.setdefault(record.provenance.get("protein_id", ""), set()).add(record.description)
        for record in records:
            if len(by_protein.get(record.provenance.get("protein_id", ""), set())) > 1:
                record.metrics["conflict"] = True
                record.provenance["conflict"] = True
