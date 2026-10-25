# PhageMine: Bacteriophage Genome Annotation and Comparative Genomics

[![CI](https://github.com/emmannaemeka/PhageMine/actions/workflows/ci.yml/badge.svg)](https://github.com/emmannaemeka/PhageMine/actions/workflows/ci.yml)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](https://bioconda.github.io/recipes/phagemine/README.html)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22794629.svg)](https://doi.org/10.5281/zenodo.22794629)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**PhageMine** is an open-source bioinformatics pipeline for bacteriophage genome
annotation, functional evidence integration and comparative genomics. It combines
phage-focused gene prediction with protein sequence and profile searches, links
functional assignments to their supporting evidence, and produces annotated
sequences, evidence tables and an HTML report.

PhageMine uses PHANOTATE for protein-coding gene prediction; PHROGs, VOGDB, Pfam
and Swiss-Prot for functional evidence; and INPHARED references for comparative
genomic context. It is designed for researchers working in phage genomics,
bacteriophage genome annotation, viral genomics and related bioinformatics
workflows.

## What it does

- Predicts protein-coding genes with PHANOTATE.
- Uses PHROGs, VOGDB, Pfam and Swiss-Prot to assess protein function.
- Reports supporting evidence, conflicting assignments and unknown functions.
- Processes individual genomes or genome collections.
- Groups related proteins into families and compares genomes with INPHARED references.
- Exports GFF3, protein FASTA and GenBank preparation files.

## Install

The current stable release is [v1.3.2](https://github.com/emmannaemeka/PhageMine/releases/tag/v1.3.2).
Python 3.10 or later is required.

Install from Bioconda, which brings the external analysis tools with it:

```bash
conda create --name phagemine --channel conda-forge --channel bioconda phagemine
conda activate phagemine
phagemine doctor
```

Or from PyPI, with the external tools installed separately:

```bash
pip install phagemine
```

To install from source instead:

```bash
git clone --branch v1.3.2 --depth 1 https://github.com/emmannaemeka/PhageMine.git
cd PhageMine
conda env create --file environment.yml
conda activate phagemine
python -m pip install .
phagemine doctor
```

Install the databases before running functional annotation:

```bash
phagemine databases install --all
```

This installs PHROGs, Swiss-Prot, Pfam and VOGDB, plus PMFDB and INPHARED
comparative references. Prepared databases require additional disk space beyond
the download sizes.

See the [installation guide](docs/installation.md) for platform requirements
and troubleshooting. The GUI is included in the normal installation:

```bash
phagemine-gui
```

Your default browser opens when the local server is ready. Keep the terminal
open while using the interface. See the [GUI guide](docs/gui.md).

### Install on Windows with WSL

You can use PhageMine on a **Windows computer through WSL**, which runs Ubuntu
underneath. The GUI opens in your normal **Windows browser**. The full pipeline
is not currently supported on native Windows; run the installation and analysis
commands inside Ubuntu, not PowerShell or Windows Command Prompt.

1. Open **PowerShell as Administrator** and install Ubuntu through WSL:

   ```powershell
   wsl --install -d Ubuntu
   ```

   Restart Windows if prompted, then open **Ubuntu** from the Start menu and
   create your Linux username and password. See Microsoft's
   [WSL installation guide](https://learn.microsoft.com/windows/wsl/install).

2. In the **Ubuntu terminal**, install Git and Miniforge. Skip the Miniforge
   download and installation if Conda is already installed inside Ubuntu:

   ```bash
   sudo apt update
   sudo apt install -y git curl ca-certificates
   curl -fL "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-$(uname -m).sh" -o /tmp/miniforge.sh
   bash /tmp/miniforge.sh -b -p "$HOME/miniforge3"
   source "$HOME/miniforge3/etc/profile.d/conda.sh"
   conda init bash
   ```

   Use Linux Conda inside Ubuntu, even if you already have Conda on Windows.
   The installer above selects your Linux CPU architecture; availability of
   PhageMine's external packages still depends on that architecture.

3. Install the current GitHub version and its dependencies **inside Ubuntu**:

   ```bash
   cd ~
   git clone https://github.com/emmannaemeka/PhageMine.git
   cd PhageMine
   conda env create --file environment.yml
   conda activate phagemine
   python -m pip install .
   phagemine doctor
   phagemine databases install --all
   phagemine doctor --deep
   ```

   Ensure the required tools and evidence resources report ready before using
   them. Database preparation needs additional disk space beyond the downloads;
   keep at least 15–20 GiB free for the complete setup.

4. Start the GUI from the same Ubuntu terminal:

   ```bash
   phagemine-gui
   ```

   If WSL does not open your Windows browser automatically, open
   [http://localhost:8501](http://localhost:8501) manually. Use the port printed
   by the launcher if it selects a different one. Keep the Ubuntu terminal open
   while using the GUI; press **Ctrl+C** to stop it.

For later sessions, open Ubuntu and run `conda activate phagemine`, followed by
`phagemine-gui`. Windows files are accessible under `/mnt/c/`; keeping large
analysis inputs and results in your Ubuntu home directory avoids cross-filesystem
I/O overhead. WSL-specific end-to-end testing is not part of the current CI matrix.

## Annotate a genome

Supply a single-genome nucleotide FASTA file:

```bash
phagemine run genome.fasta --output results/genome --threads 8
```

Open `results/genome/report.html` to browse the results. The output directory
also contains `annotation.tsv`, `genes.gff3`, `proteins.faa` and supporting
evidence files.

To check the installation with the bundled example:

```bash
phagemine run examples/demo_phage.fasta --output results/example
```

Installed databases are detected automatically. Without them, the example
checks gene prediction and report generation; functional evidence is unavailable.

Run `phagemine --help` for all commands and `phagemine run --help` for annotation
options.

## Reading the results

Product names describe the function supported by the available evidence.
Proteins without a supported function remain hypothetical or are identified as
conserved proteins of unknown function. Gene symbols, protein products and
locus identifiers are reported separately.

PhageMine is research software. Functional predictions require review,
and confidence labels describe evidence strength rather than calibrated
probabilities. See the [limitations](docs/limitations.md) for interpretation.

## Citation and support

Citation information is available in [CITATION.cff](CITATION.cff) and the
[archived release](https://doi.org/10.5281/zenodo.22794629).

Report bugs or request features through
[GitHub Issues](https://github.com/emmannaemeka/PhageMine/issues).
PhageMine is distributed under the [MIT licence](LICENSE).
