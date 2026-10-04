# PhageMine

PhageMine annotates bacteriophage genomes and organizes the evidence behind
protein function predictions. It combines phage gene prediction with protein
sequence and profile searches, and produces annotated sequences, tables and
an HTML report.

## What it does

- Predicts protein-coding genes with PHANOTATE.
- Uses PHROGs, VOGDB, Pfam and Swiss-Prot to assess protein function.
- Reports supporting evidence, conflicting assignments and unknown functions.
- Processes individual genomes or genome collections.
- Groups related proteins into families and compares genomes with INPHARED references.
- Exports GFF3, protein FASTA and GenBank preparation files.

## Install

Python 3.10 or later is required. Use Conda to install the external analysis tools:

```bash
git clone https://github.com/emmannaemeka/PhageMine.git
cd PhageMine
conda env create --file environment.yml
conda activate phagemine
python -m pip install .
phagemine doctor
```

Install the annotation databases before running functional annotation:

```bash
phagemine databases install pfam
phagemine databases install vogdb
phagemine databases install swissprot
phagemine databases install phrogs
```

For protein-family and reference-genome comparisons, also install `pmfdb` and
`inphared`. Database preparation requires additional disk space; allow at least
15–20 GiB for the complete resource set.

See the [installation guide](docs/installation.md) for platform requirements
and troubleshooting.

## Annotate a genome

Supply a single-genome nucleotide FASTA file:

```bash
phagemine run genome.fasta --output results/genome --evidence full --threads 8
```

Open `results/genome/report.html` to browse the results. The output directory
also contains `annotation.tsv`, `genes.gff3`, `proteins.faa` and supporting
evidence files.

To check the installation with the bundled example:

```bash
phagemine run examples/demo_phage.fasta --output results/example
```

This example runs without evidence databases and checks gene prediction and
report generation. Use `--evidence full` with the installed databases for
functional annotation.

Run `phagemine --help` for all commands and `phagemine run --help` for annotation
options.

## Reading the results

Product names describe the function supported by the available evidence.
Proteins without a supported function remain hypothetical or are identified as
conserved proteins of unknown function. Gene symbols, protein products and
locus identifiers are reported separately.

PhageMine is beta research software. Functional predictions require review,
and confidence labels describe evidence strength rather than calibrated
probabilities. See the [limitations](docs/limitations.md) for interpretation.

## Citation and support

Citation information is available in [CITATION.cff](CITATION.cff) and the
[archived release](https://doi.org/10.5281/zenodo.22794629).

Report bugs or request features through
[GitHub Issues](https://github.com/emmannaemeka/PhageMine/issues).
PhageMine is distributed under the [MIT licence](LICENSE).
