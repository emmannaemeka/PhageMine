# Installation

Install PhageMine in a clean environment. A Conda installation is recommended
because it supplies PHANOTATE, Prodigal, HMMER, MMseqs2, DIAMOND and Mash together:

```bash
conda create -n phagemine phagemine \
  --channel conda-forge \
  --channel bioconda \
  --strict-channel-priority
conda activate phagemine
```

The Bioconda recipe for PhageMine 1.0.2 requires PHANOTATE, Prodigal, HMMER,
MMseqs2, DIAMOND and Mash. Do not treat package installation alone as proof that these
compiled programs run on the host. `phagemine doctor` executes each program
and reports `BROKEN` when a version probe fails, including dynamic-linker
errors.

For a local wheel:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install phagemine-1.0.0-py3-none-any.whl
```

## Required post-installation database setup

The software package does not bundle the large biological evidence databases.
Immediately after installing PhageMine, run:

```bash
phagemine databases install --all
phagemine doctor
```

All required executables and requested capabilities must report `READY` before
a production full-evidence run. `BLOCKED/INVALID` means that a database is
registered but cannot currently operate, commonly because a required
executable is missing or broken; it does not mean the database must be
downloaded again.

The Conda recipe uses the supplied `packaging/bioconda/post-link.sh` reminder.
The same instructions appear when `phagemine` or `phagemine --help` is run.
`phagemine doctor` repeats actionable installation commands for every missing
resource.

Install resources separately when bandwidth or storage is constrained:

```bash
phagemine databases install pfam
phagemine databases install vogdb
phagemine databases install swissprot
phagemine databases install phrogs
phagemine databases install pmfdb
phagemine databases install inphared
phagemine doctor
```

Production gene prediction requires the external PHANOTATE script. Pass its
path explicitly when it is not available as `phanotate` on `PATH`, for example
`--phanotate /path/to/phanotate.py`.

Minimal Core run:

```bash
phagemine run genome.fasta \
  --phanotate /path/to/phanotate.py
```

Without `--output`, PhageMine writes
`./genome_phagemine_results` in the directory where the command is run. Use an
absolute `--output` path when a different location is required. The complete
`run` workflow also compares PHANOTATE and Prodigal calls and writes review
tables; it does not automatically discard discordant ORFs.

Evidence databases are not required for Core. Missing resources are reported
as unavailable and the Core workflow remains usable, but full evidence mode
requires the four annotation resources to be `READY` in `phagemine doctor`.
PMFDB and INPHARED enable the two comparative discovery capabilities. The
complete six-resource download is approximately 3.4 GiB; keep at least
15–20 GiB free during preparation.
