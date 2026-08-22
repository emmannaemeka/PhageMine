# External databases and executables

PhageMine Core does not bundle Pfam, VOGDB, PHROGs, Swiss-Prot, or PMFDB.
These independently installed, versioned resources are supplied only when
evidence-enhanced annotation or external PMF validation is requested.

Optional components include HMMER (`hmmscan`) with Pfam/VOGDB HMM files,
PHROGs profiles and annotations with MMseqs2, DIAMOND with Swiss-Prot and
matching metadata, and MMseqs2 with a PMFDB target database/index.

Record paths, versions, and checksums in run provenance. Missing optional
resources are recorded as unavailable, not as biological absence. PMFDB is an
independently versioned external knowledgebase.

## Automated installation

```bash
phagemine databases install pfam
phagemine databases install vogdb
phagemine databases install swissprot
phagemine databases install phrogs

# Equivalent complete installation
phagemine databases install --all
phagemine doctor
```

The default database root is the platform user-data directory:

- macOS: `~/Library/Application Support/PhageMine/databases`
- Linux: `${XDG_DATA_HOME:-~/.local/share}/phagemine/databases`
- Windows: `%LOCALAPPDATA%\PhageMine\databases`

Use `--directory PATH` to select another disk. `--dry-run` displays the
estimated compressed download before network access. The complete download is
approximately 2.4 GiB; preparation requires more temporary and final storage.

## Installation behavior and provenance

| Resource | Distribution | Preparation | Registered sidecar |
|---|---|---|---|
| Pfam | Provider-described current Pfam release | decompress and `hmmpress` | release and checksum files |
| VOGDB | Pinned VOGDB release 235 | deterministic HMM concatenation and `hmmpress` | `vog.annotations.tsv.gz` |
| Swiss-Prot | Provider-described current UniProtKB release | `diamond makedb` | `uniprot_sprot.dat.gz` |
| PHROGs | PHROGs v4 from the checksum-pinned Pharokka 1.11.0 distribution | extract MMseqs2 profile database | `phrog_annot_v4.tsv` |

Each completed directory contains `install_manifest.json` with release,
provider, source URL, checksum and preparation information. Registration occurs
only after preparation. A resource is reported as installed only after the
normal PhageMine validator returns `READY`.

Downloads use a `.part` file and HTTP range requests when the provider supports
them, allowing an interrupted download to resume. Successfully used archives
are removed unless `--keep-downloads` is supplied.

Existing locally prepared snapshots remain supported through `phagemine
databases register`. Automatic updates are intentionally not silent: rerun an
installer with `--force` only after reviewing the new release and preserving
the older manifest needed for reproducibility.
