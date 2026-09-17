# Installation

## Quick start

Use the repository environment for a first installation. It installs the
external bioinformatics executables from Conda channels and leaves the Python
package installation explicit:

```bash
git clone https://github.com/emmannaemeka/PhageMine.git
cd PhageMine
conda env create --file environment.yml
conda activate phagemine
python -m pip install .
phagemine --help
phagemine doctor
phagemine databases install --all
phagemine doctor --deep
```

The command sequence assumes a supported Conda-compatible client. The package
and unit-test CI matrix covers Ubuntu Linux, Intel macOS (`macos-13`) and
Apple-silicon macOS (`macos-14`); external-tool integration is exercised on
Ubuntu. This is not a claim that every Conda package is available on every
local architecture.

The environment supplies PHANOTATE, Prodigal, HMMER (`hmmscan`), MMseqs2,
DIAMOND, Mash and BLASTN. PyHMMER is a Python library used by the PHROGs
provider; it is not the HMMER executable. `table2asn` is optional. Database
downloads do not install missing executables.

For development, install the test extras in editable mode after the environment
has been created:

```bash
python -m pip install -e ".[test]"
```

That developer-only extra also installs the pinned Pyrodigal, Prodigal-gv and
Pyrodigal-rv Python providers used by diagnostic, RNA and test workflows. They
are optional; the default PHANOTATE-based Core workflow does not require them.

If `pip` reports that neither `setup.py` nor `pyproject.toml` exists, the
command is not running inside `PhageMine`; run `cd /path/to/PhageMine` first.
If `phagemine` is not found, activate the environment with `conda activate
phagemine` and rerun `python -m pip install .`. A new Terminal session also
requires `cd /path/to/PhageMine` and `conda activate phagemine` again.

The Bioconda recipe for PhageMine 1.2.0 must require PHANOTATE, Prodigal, HMMER,
MMseqs2, PyHMMER, DIAMOND, Mash and BLASTN (the Bioconda `blast` package). Do not treat package installation alone as proof that these
compiled programs run on the host. `phagemine doctor` executes each program
and reports `BROKEN` when a version probe fails, including dynamic-linker
errors.

## Required post-installation database setup

The software package does not bundle the large biological evidence databases.
Immediately after installing PhageMine, run:

```bash
phagemine databases install --all
phagemine doctor
```

For release validation or diagnosis of a database that appears ready but fails
at runtime, use `phagemine doctor --deep`. It opens the registered PHROGs
MMseqs2 and PyHMMER databases and reports format errors before an analysis run.

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
INPHARED workflow uses Mash only to shortlist references and BLASTN to compute
PhageMine's bidirectional length-normalized nucleotide similarity; this is not
presented as VIRIDIC output. Mash distance is
never converted into percentage similarity. The
complete six-resource download is approximately 3.4 GiB; keep at least
15–20 GiB free during preparation.
