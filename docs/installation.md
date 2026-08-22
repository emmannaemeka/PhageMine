# Installation

Install PhageMine in a clean environment. A Conda installation is recommended
because it supplies PHANOTATE, HMMER, MMseqs2, DIAMOND and Mash together:

```bash
conda create -n phagemine phagemine \
  --channel conda-forge \
  --channel bioconda \
  --strict-channel-priority
conda activate phagemine
```

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
phagemine run genome.fasta --output results/run \
  --phanotate /path/to/phanotate.py
```

Evidence databases are not required for Core. Missing resources are reported
as unavailable and the Core workflow remains usable, but full evidence mode
requires the four annotation resources to be `READY` in `phagemine doctor`.
PMFDB and INPHARED enable the two comparative discovery capabilities. The
complete six-resource download is approximately 3.4 GiB; keep at least
15–20 GiB free during preparation.
