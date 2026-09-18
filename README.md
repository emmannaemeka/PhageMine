# PhageMine

[![Release: v1.2.0](https://img.shields.io/badge/release-v1.2.0-2ea44f)](https://github.com/emmannaemeka/PhageMine/releases/tag/v1.2.0)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22794629.svg)](https://doi.org/10.5281/zenodo.22794629)

PhageMine v1.2 is a phage-focused genome annotation workflow that combines
PHANOTATE structural gene prediction with evidence-supported functional
annotation. It keeps gene calls, evidence, uncertainty, provenance, and
review flags separate so that a named product is never mistaken for
experimental confirmation.

This repository contains the source, tests, installation metadata, examples,
and compact benchmark summaries suitable for public reproducibility. Large
evidence databases and local benchmark evidence are intentionally kept out of
the public source tree.

## Local graphical interface (GUI v0.1)

Install the optional GUI dependencies and launch the interface:

```bash
pip install "phagemine[gui]"
phagemine-gui
```

`phagemine gui` is also supported. The Streamlit interface runs locally, does
not upload analysis files to an external server, and calls the same PhageMine
analysis engine used by the CLI. The CLI remains fully supported and its core
installation does not depend on Streamlit.

GUI documentation assets, including a future real screenshot, belong in
`docs/images/phagemine-gui/`. No screenshot is included until one has been
captured from a validated release build.

PhageMine is beta research software for evidence-based bacteriophage genome annotation and discovery mining. It annotates what can be supported by evidence and organizes what remains unknown across a cohort. Unknown does not mean novel, and its rule-based evidence-strength labels are not calibrated probabilities.

## Complete workflow

```text
Genome FASTA
    ↓
Phage-focused CDS prediction
    ↓
Functional annotation
PHROGs • VOGDB • Pfam • Swiss-Prot
    ↓
Evidence-supported annotated genome
    ↓
INPHARED comparative analysis
    ↓
Nearest known phages + genomic/taxonomic context
```

INPHARED adds comparative genomic context through Mash nearest-reference screening
and, when configured, optional bidirectional BLASTN confirmation. It reports
supported reference and host/taxonomy metadata; it is not formal ICTV
classification.

## Why PhageMine?

Phage genomes contain many hypothetical or uncharacterized proteins. Conventional annotation often stops at “hypothetical protein”. PhageMine combines conservative evidence fusion for defensible annotation with cohort-level discovery of recurrent, context-preserved protein families. Predictions are computational hypotheses, not experimental confirmation.

## Key features

- PHANOTATE CDS prediction and genome validation
- Pfam, VOGDB, Swiss-Prot, and dual-backend PHROGs evidence (MMseqs2 plus PyHMMER)
- deterministic evidence fusion, functional state, proposed function, and confidence
- genomic context, modules, candidate ranking, QC, and GenBank pre-submission files
- batch annotation and HTML reports with linked protein records
- PMF (PhageMine protein family) clustering, recurrence, synteny/context analysis, and PMFDB validation
- publication-oriented PNG/SVG figures and retained figure source tables
- checkpointing, fingerprints, resume, provenance, and BOTH-mode evidence reuse

## Benchmark Results — v1.2

PhageMine v1.2 was evaluated on seven curated reference bacteriophage genomes
against Pharokka and Prokka.

| Tool | Predicted CDSs | Strict F1 | Relaxed F1 | Named products |
|------|---------------:|----------:|-----------:|---------------:|
| PhageMine | 804 | 0.7762 | 0.8821 | 439 |
| Pharokka | 804 | 0.7762 | 0.8821 | 398 |
| Prokka | 675 | 0.8846 | 0.9352 | 328 |

PhageMine and Pharokka produced identical CDS coordinates across the
seven-genome benchmark panel. Of PhageMine's 216 strict non-exact predictions,
141 were alternative-boundary or same-strand-overlap cases. Consequently,
PhageMine's F1 increased from 0.7762 under strict exact-coordinate scoring to
0.8821 under relaxed gene-level scoring. Prokka showed stronger
reference-relative structural agreement on this curated seven-phage panel.
PhageMine assigned 439 named products compared with 398 for Pharokka and 328
for Prokka. Named-product yield is not equivalent to functional annotation
accuracy. These results apply to this seven-genome curated benchmark and must
not be interpreted as universal tool rankings. INPHARED comparative/taxonomic
performance was not evaluated by this structural benchmark.

![Strict versus relaxed F1 across tools](docs/benchmark_v1.2/figures/strict_vs_relaxed_f1.png)

*Strict and relaxed reference-relative structural F1. See the full benchmark
page for per-genome results and definitions.*

![Named-product yield comparison](docs/benchmark_v1.2/figures/functional_yield.png)

*Named-product yield is an output metric and is not functional accuracy. The
previous erroneous per-genome PhageMine zeros have been corrected from frozen
annotation tables; the aggregate remains 439. See the [correction note](docs/BENCHMARK.md)
and [functional reference-concordance analysis](docs/FUNCTIONAL_BENCHMARK.md).*

**Full benchmark methodology, per-genome results, figures and limitations: [docs/BENCHMARK.md](docs/BENCHMARK.md)**

## The three modes

### Annotation mode

`genome FASTA → validation → PHANOTATE → proteins → Pfam/VOGDB/Swiss-Prot/PHROGs → evidence fusion → classification → annotated genome/proteins → INPHARED comparative context → ranking/QC → figures/report → GenBank package`

For DNA genomes, PHANOTATE is the sole primary structural caller. Pyrodigal and
Prodigal-gv are diagnostic corroboration/conflict-detection callers; they never
vote, alter PHANOTATE boundaries, or add caller-specific CDSs. Evidence
annotates predicted proteins and never silently changes primary ORF boundaries.
A report records coordinates, strand, classification, proposed function,
confidence, supporting evidence, and a deterministic reason.

Illustrative example (not a guaranteed result):

```text
Protein: PM_000023
Genome: NC_011107.1
Classification: PROBABLE_FUNCTION
Proposed function: Terminase large subunit
Evidence: Pfam terminase ATPase; VOGDB terminase; PHROGs DNA packaging;
          Swiss-Prot terminase homolog
Confidence: HIGH
Reason: Multiple independent evidence sources support a DNA-packaging
        terminase assignment.
```

### Discovery mode

`genomes → cohort QC → PHANOTATE per genome → protein pooling → exact AA deduplication → pooled evidence → remapping → classification → PMF clustering → recurrence → genomic context/synteny → PMFDB → ranking → figures/report`

A PMF is a homologous PhageMine protein family based on biological sequence similarity. Exact amino-acid deduplication is only a computational optimization and does not define a PMF. PMFDB is the PhageMine protein-family knowledgebase. Its conservative results distinguish experimentally characterized, predicted-function, uncharacterized, mixed, and metadata-incomplete matches. `NO_EXTERNAL_MATCH` does not mean novel, and an INPHARED product label is retained as a computational prediction rather than experimental evidence.

### Both mode

Both mode completes per-genome annotation, then reuses valid proteins, coordinates, evidence, classifications, and provenance for cohort discovery. Equivalent expensive searches are not intentionally repeated.

## Quick start

PhageMine supports Python 3.10 or newer. The recommended first-time route is
Conda (or a compatible `mamba`/`micromamba` client), because the repository's
`environment.yml` installs the external command-line tools alongside Python.
Run these commands exactly from a Terminal:

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

The first two Doctor commands check the executable environment before database
downloads. The final command checks both executable availability and the
operational database formats. The complete database download is approximately
3.4 GiB compressed; preparation and temporary files require additional space,
so keep at least 15–20 GiB free. Database installation is a separate step: it
does not install missing executables.

The Conda environment installs these external executables: PHANOTATE,
Prodigal, HMMER (`hmmscan`), MMseqs2 (`mmseqs`), DIAMOND (`diamond`), Mash
(`mash`) and BLASTN (`blastn`). PyHMMER is a Python library used by the PHROGs
provider, not an executable replacement for HMMER. `table2asn` is optional and
is not required for ordinary analysis. Mash selects candidate references;
BLASTN supplies PhageMine's bidirectional, length-normalized nucleotide
comparison. It is not presented as VIRIDIC output and does not assign taxonomy.

The repository's automated matrix covers Ubuntu Linux, Intel macOS (`macos-15-intel`)
and Apple-silicon macOS (`macos-14`) for the Python package and unit suite. The
external-tool Conda workflow is exercised on Ubuntu. Local Conda solver and
Bioconda availability can vary by architecture; if the environment solver
cannot provide a package, use the documented error and platform-specific
package channel rather than a universal command that has not been tested.

For developers who need source changes to take effect immediately, use the
separate editable route after creating the environment:

```bash
python -m pip install -e ".[test]"
```

Do not use editable installation for a normal user installation.

The native Pyrodigal, Prodigal-gv and Pyrodigal-rv Python providers are
developer/test extras. They provide optional diagnostic and RNA-provider
capabilities; they are not needed for the default PHANOTATE-based Core run.
The developer command above installs their pinned versions for the test suite.

### Troubleshooting installation

* `Neither setup.py nor pyproject.toml found` means the command was run outside
  the project directory. Run `cd PhageMine` (or `cd /path/to/PhageMine`) and
  repeat `python -m pip install .`.
* `phagemine: command not found` usually means the environment is not active.
  Run `conda activate phagemine`; then check `python -m pip show phagemine`.
  If it is still absent, repeat `python -m pip install .` while the environment
  is active.
* A Doctor message about a missing external tool means the executable is not
  available or is broken. It is not fixed by downloading a database. Activate
  `phagemine`, verify `command -v hmmscan mmseqs diamond mash blastn phanotate`
  (and `prodigal`), then run `phagemine doctor` again.
* A new Terminal session does not retain activation. Run `cd /path/to/PhageMine`
  followed by `conda activate phagemine` before using `phagemine`.

For an operational database-format check before a long run, use:

```bash
phagemine doctor --deep
```

This opens the registered PHROGs MMseqs2 and PyHMMER databases and fails if a
database is only superficially present but cannot actually be read.

Do not proceed from executable presence alone. Doctor requires each executable
probe to exit successfully and reports runtime/linker failures as `BROKEN`.
An installed database that is unusable because its required executable is
missing or broken is shown as `BLOCKED/INVALID`, not as an instruction to
download the database again.

The database installer reports an approximately 3.4 GiB compressed download
before it begins. Prepared databases and temporary working files require
additional disk space; keep at least 15–20 GiB free for `--all`.

## Database configuration

Large databases are not bundled. Install one database at a time when preferred:

```bash
phagemine databases install pfam
phagemine databases install vogdb
phagemine databases install swissprot
phagemine databases install phrogs
phagemine databases install pmfdb
phagemine databases install inphared
phagemine doctor
```

Each installer downloads a provider release, verifies the published checksum
where available, prepares the searchable database, writes an installation
manifest, registers the resource, and requires it to pass validation before
reporting `READY`. Downloads can be resumed from the persistent `.downloads`
directory. Use `--directory`, `--force`, `--keep-downloads`, `--dry-run`, or
`--json` as needed.

Researchers with an existing local snapshot can register it without downloading:

```bash
phagemine databases register pfam /path/to/Pfam-A.hmm --name Pfam-A
phagemine databases register vogdb /path/to/VOGDB.hmm --name VOGDB \
  --annotations /path/to/vog.annotations.tsv.gz
phagemine databases register swissprot /path/to/uniprot_sprot.dmnd \
  --name Swiss-Prot --metadata /path/to/uniprot_sprot.dat.gz
phagemine databases register phrogs /path/to/phrogs_profile_db \
  --name PHROGs --annotations /path/to/phrog_annotations.tsv
phagemine databases register pmfdb /path/to/PMFDB-INPHARED-2026-04-07 \
  --name PMFDB-INPHARED
phagemine databases register inphared_genomes /path/to/reference_phage_genomes.fna \
  --name INPHARED-Genomes --metadata /path/to/genome_metadata.tsv
```

`pmfdb` converts the pinned 7 April 2026 INPHARED proteins, protein-to-genome
mapping and host/taxonomy table into a versioned PMFDB release and builds its
MMseqs2 index. `inphared` prepares the corresponding genome FASTA and Mash
sketch. INPHARED product names remain computational predictions and are never
promoted to experimental characterization. Check readiness with:

```bash
phagemine databases status
phagemine doctor --json
```

`READY` means the registered path and required sidecars/tools are available. An unavailable resource is reported conservatively and its evidence is not fabricated.

## Quick start

Install the pinned INPHARED comparative resource, check the environment, and
run one genome:

```bash
phagemine databases install inphared
phagemine doctor
phagemine run genome.fasta --output results/genome
```

The run writes predicted CDSs and proteins, functional product assignments,
evidence and provenance records, INPHARED nearest-phage results when the
resource is available, comparative reports, and a run manifest. See
[docs/INPHARED.md](docs/INPHARED.md) for comparative outputs and
[docs/BENCHMARK.md](docs/BENCHMARK.md) for validated benchmark interpretation.

```bash
# Optional: annotate a batch with a different destination.
phagemine batch genomes/ --output results/annotation --mode annotate \
  --evidence full --threads 8

phagemine batch genomes/ --output results/discovery --mode discover \
  --evidence full --threads 8

phagemine batch genomes/ --output results/both --mode both \
  --evidence full --threads 8
```

## Interpreting annotation

The main table uses plain-language classifications such as “Specific function strongly supported”, “Likely function supported by evidence”, “Protein domain detected; full function unknown”, “Conserved in phages; function unknown”, and “No reliable function identified”. Proposed function answers “what does PhageMine think it does?”, while rule-based evidence strength and best evidence explain why. These strength categories are explicitly uncalibrated, and domain notes are separated from product names. Machine-readable states remain in `functional_classification.tsv`. A prediction is not experimental confirmation. PhageMine does not infer lifestyle.

Batch tables include protein ID, coordinates, strand, classification, proposed function, confidence, evidence summaries, and PMF where available. Protein IDs link to detailed records rather than reducing a row to a state label alone.

## Protein retrieval

```bash
phagemine extract protein PM_000023 --run results/genome --protein-fasta
phagemine extract protein PM_000023 --run results/genome --cds-fasta
phagemine extract protein PM_000023 --run results/genome --evidence
```

Each record also links genomic context and, in discovery results, PMF membership. The internal PhageMine protein ID is the stable key connecting the report row, amino-acid FASTA, nucleotide CDS, evidence, and context.

## Important outputs

Annotation outputs include the concise `annotation.tsv`, detailed `functional_classification.tsv/json`, `scientific_validation_status.json`, `evidence.json`, `candidate_ranking.tsv`, ID-only `proteins.faa`, product-labelled `annotated_proteins.faa`, `cds.fna`, `genes.gff3`, `gene_call_confidence.tsv`, `gene_calls_for_review.tsv`, `genomic_context.tsv/json`, `modules.tsv/json`, `quality_control.json`, `report.html`, `report.md`, `protein_details/`, `figures/`, `figure_data/`, and `genbank_submission/`.

`hallmark_completeness.tsv` is a backward-compatible annotation-text screen for
major capsid, portal, terminase, tail, tape-measure, replication and lysis
components. It is not an independent HMM search or genome-completeness
estimate. `NOT_ESTABLISHED` never means biological absence.
`annotation_review.tsv` is the short manual-curation queue combining uncertain
functions, disputed gene models and unresolved hallmark components. Detailed
classification records also state the evidence tier used for each conclusion.

Discovery outputs include `pmf_families.tsv`, `pmf_members.tsv`, `family_recurrence.tsv`, `cohort_unknown_proteome.tsv`, `conserved_neighbourhoods.tsv`, `pmfdb_validation.tsv`, `discovery_ranking.tsv`, `discovery_report.html`, `discovery_report.md`, `figures/`, `figure_data/`, and `pooled_manifest.json`.

When the INPHARED genome resource is READY, discovery also writes
`inphared_nearest_phages.tsv` as the accession-level audit trail and
`inphared_summary.tsv` as the duplicate-collapsed researcher summary.

## What a researcher receives from `phagemine run`

Run a genome with:

```bash
phagemine run genome.fasta --output results/genome
```

The workflow validates the input, calls CDSs with PHANOTATE, translates the
CDSs, searches the registered evidence resources, and writes a linked set of
machine-readable and human-readable outputs. The primary downstream files are
`genes.gff3` (coordinates and strand), `proteins.faa` (stable protein IDs),
`cds.fna` (nucleotide CDSs), `annotation.tsv` (one row per predicted CDS),
`functional_classification.tsv/json` (structured evidence and classifications),
`evidence.json` (source-level evidence), `gene_call_confidence.tsv`,
`quality_control.json`, and `report.html`/`report.md`.

An annotation row contains the protein ID, start, end, strand, amino-acid
length, gene label when available, product or proposed function, confidence,
classification, evidence sources, and review flags. A small illustrative row is:

```text
protein_id  start  end   strand  product                    confidence  evidence_sources
PM_000023   1452   2387  +       terminase large subunit     HIGH        PHROGs;VOGDB;Swiss-Prot
```

Coordinates identify a computational CDS hypothesis. A product assignment is
an evidence-backed annotation statement, not a laboratory result. “Hypothetical”
and “uncharacterized” remain meaningful outputs when evidence does not support
greater specificity.

PHANOTATE is the primary DNA caller. The supported small-genome rule selects
Prodigal metagenomic mode for sequences below its single-genome minimum when
secondary corroboration is requested; secondary callers do not silently alter
the PHANOTATE primary coordinates. Functional evidence may come from PHROGs,
VOGDB, Pfam, Swiss-Prot, and related registered resources. Each output records
which sources contributed to the conclusion.

## System requirements and troubleshooting

Use Python 3.10 or newer. A Conda environment is recommended for PHANOTATE,
Prodigal, HMMER, MMseqs2, DIAMOND, Mash, BLASTN, PyHMMER, and the optional
Streamlit GUI. Core annotation can run without evidence databases, but full
functional annotation requires the relevant resources to be installed and
reported `READY` by `phagemine doctor`.

If a run stops before annotation, run `phagemine doctor` and then
`phagemine doctor --deep`. Check executable versions, database manifests,
write permissions, temporary disk space, and the run's `scientific_validation_status.json`.
Keep 15–20 GiB free when installing the full evidence set. Do not interpret an
unavailable database as zero biological evidence.

## Benchmark summary

The corrected seven-reference-phage benchmark and its limitations are
documented in [docs/BENCHMARK.md](docs/BENCHMARK.md). It reports measured
strict and relaxed structural metrics, functional annotation yield, and the
distinction between the frozen original benchmark and the corrected analysis.
The benchmark does not establish universal superiority for any tool.

## INPHARED comparative genomic context

When a validated INPHARED resource is installed, `phagemine run` and
`phagemine annotate` append nearest-reference genomic context automatically;
batch `discover` and `both` workflows do the same for their discovery output.
Install it with:

```bash
phagemine databases install inphared
phagemine doctor
phagemine run genome.fasta --output results/genome
```

The pinned INPHARED release is **2026-04-07**. Installation verifies provider
checksums, prepares the reference FASTA and metadata sidecars, and builds a
per-genome Mash sketch. Mash ranks nearest reference phages by distance,
p-value, and matching hashes. When `blastn` and reference sequences are
available, PhageMine also reports its own bidirectional, length-normalized
nucleotide similarity and alignment coverage. This is not ANI or VIRIDIC.
Reference accession, description, host genus, phage genus, subfamily, and
family are returned as metadata; PhageMine does not make formal ICTV
assignments. Partial or poorly aligned queries are withheld from numerical
taxonomy, and a 95% species working boundary is labelled as requiring ICTV
review.

Results are written to `comparative/inphared_nearest_phages.tsv` and `.json`,
`comparative/inphared_summary.tsv`, the run manifest, and the HTML/Markdown
report. If the resource is unavailable or invalid, the workflow records
`SKIPPED` with a reason and does not fabricate zero matches. See
[docs/INPHARED.md](docs/INPHARED.md) for columns, provenance, interpretation,
and limitations. The structural benchmark in [docs/BENCHMARK.md](docs/BENCHMARK.md)
did not validate INPHARED taxonomic performance.

## Limitations

The validation panel contains seven curated phages and does not represent the
full diversity of phage genomes. Exact-coordinate scoring penalizes legitimate
alternative starts and stops; relaxed overlap scoring addresses gene-level
detection but does not prove biological truth. Named-product yield is not
functional accuracy. Runtime ranking is unsupported by the available timing
provenance. The historical T4 PHANOTATE 294-versus-297 discrepancy remains
unresolved. Independent validation on broader panels is appropriate.

## Reproducibility and citation

Pin the release version, record `phagemine doctor --json`, retain database
manifests and checksums, and preserve the complete output directory. The
compact public benchmark package includes accession lists, reference
provenance, tables, figures, and checksums; third-party databases are not
redistributed. PhageMine v1.2.0 is permanently archived on Zenodo.

### Citation

Nnadi, Nnaemeka Emmanuel. (2026). *PhageMine v1.2.0*. Zenodo.
https://doi.org/10.5281/zenodo.22794629

See [CITATION.cff](CITATION.cff) for machine-readable citation metadata and
[CHANGELOG.md](CHANGELOG.md) for release history. Contributions and issue
reports are welcome through the GitHub repository.
For zero-distance Mash hits, PhageMine directly compares the query and
reference nucleotide sequences, including reverse-complement and
rotation-equivalent representations. Numerical-boundary interpretation is
withheld for partial or poorly aligned queries. Formal VIRIDIC/phylogenetic
analysis and current family-specific ICTV criteria remain necessary; taxonomy
is not inferred.

## Independent annotation comparison

PhageMine can import public standard-tool outputs without running or modifying
those tools:

```bash
phagemine benchmark --output comparison \
  --phagemine-results ./genome_phagemine_results \
  --pharokka-gff /path/to/pharokka.gff \
  --phold-genbank /path/to/phold.gbk \
  --multiphate-gff /path/to/multiPhATE2.gff \
  --prokka-gff /path/to/prokka.gff

# Accuracy metrics require an expert-reviewed truth set:
phagemine benchmark --output validated-comparison \
  --phagemine-results ./genome_phagemine_results \
  --truth-genbank /path/to/expert_reviewed_truth.gbk
```

The comparison reports coordinate agreement separately from product-name
agreement. Without `--truth-gff` or `--truth-genbank`, the manifest is labelled
`TOOL_AGREEMENT_ONLY` and cannot support accuracy or superiority claims. See
`docs/validation-protocol.md`.
It also writes `hallmark_comparison.tsv` and a checksummed
`benchmark_manifest.json`. No structural-search software is required.

The TSV files are stable, machine-readable tables; JSON files preserve detailed evidence/provenance; reports provide navigable summaries and downloads.

## Figures

Annotation figures include genome functional maps, functional-state summaries, evidence-support matrices, and functional-category distributions where the data support them. Batch output includes genome-by-functional-state comparison. Discovery output includes functional distribution, PMF recurrence heatmap/distribution, discovery ranking, PMFDB source data, and high-priority genomic-context diagrams. Plotting source tables are retained in `figure_data/` for reproducibility.

## Resume and checkpoints

Runs store stage checkpoints, input fingerprints, checksums, database provenance, parameters, and tool versions. Resume reuses valid completed stages and reruns only missing or stale stages plus required downstream dependencies. In BOTH mode, valid annotation evidence is reused for discovery; incompatible evidence is never silently reused.

## GenBank pre-submission

The `genbank_submission/` package contains genome FASTA, CDS/protein outputs, feature tables, supplied metadata/provenance, and validation results. PhageMine does not invent missing metadata. `table2asn` can be used when installed, but is optional.

## Conservative scientific interpretation

- Unknown does not mean novel.
- Unresolved does not mean biologically unimportant.
- `NO_EXTERNAL_MATCH` does not prove novelty.
- Genomic context is supporting evidence, not proof of function.
- Sequence similarity supports homology but does not automatically establish biochemical function.
- Experimental validation remains necessary for biological claims.

## Development and citation

Run the test suite with `PYTHONPATH=src pytest -q`. See
`docs/release-checklist.md` before publishing an artifact, and see
`CHANGELOG.md`, `CITATION.cff`, and `LICENSE` for release metadata and terms.
