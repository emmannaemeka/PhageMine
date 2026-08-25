# PhageMine

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

PhageMine is an evidence-based bacteriophage genome annotation and discovery-mining platform. It annotates what can be supported by evidence and organizes what remains unknown across a cohort. Unknown does not mean novel.

## Why PhageMine?

Phage genomes contain many hypothetical or uncharacterized proteins. Conventional annotation often stops at “hypothetical protein”. PhageMine combines conservative evidence fusion for defensible annotation with cohort-level discovery of recurrent, context-preserved protein families. Predictions are computational hypotheses, not experimental confirmation.

## Key features

- PHANOTATE CDS prediction and genome validation
- Pfam, VOGDB, Swiss-Prot, and PHROGs evidence
- deterministic evidence fusion, functional state, proposed function, and confidence
- genomic context, modules, candidate ranking, QC, and GenBank pre-submission files
- batch annotation and HTML reports with linked protein records
- PMF (PhageMine protein family) clustering, recurrence, synteny/context analysis, and PMFDB validation
- publication-oriented PNG/SVG figures and retained figure source tables
- checkpointing, fingerprints, resume, provenance, and BOTH-mode evidence reuse

## The three modes

### Annotation mode

`genome FASTA → validation → PHANOTATE → proteins → Pfam/VOGDB/Swiss-Prot/PHROGs → evidence fusion → classification → context/modules → ranking/QC → figures/report → GenBank package`

Evidence annotates predicted proteins; it never silently changes PHANOTATE ORF boundaries. A report records coordinates, strand, classification, proposed function, confidence, supporting evidence, and a deterministic reason.

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

## Installation

PhageMine supports Python 3.10 or newer.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

Production runs require these executables on `PATH`: PHANOTATE, HMMER
(`hmmscan`), MMseqs2 (`mmseqs`), DIAMOND (`diamond`), Mash (`mash`) and BLASTN
(`blastn`) when whole-genome INPHARED comparison is installed. Mash selects
candidate references; BLASTN supplies VIRIDIC-compatible intergenomic
similarity. `table2asn` is optional and
never blocks the normal GenBank pre-submission package.

After installation, PhageMine displays the database setup commands in its
top-level help. Install all evidence databases and validate the environment:

```bash
phagemine databases install --all
phagemine doctor
```

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

```bash
phagemine run genome.fasta

# Optional: choose a different destination explicitly.
phagemine run genome.fasta --output /path/to/results/genome

phagemine batch genomes/ --output results/annotation --mode annotate \
  --evidence full --threads 8

phagemine batch genomes/ --output results/discovery --mode discover \
  --evidence full --threads 8

phagemine batch genomes/ --output results/both --mode both \
  --evidence full --threads 8
```

## Interpreting annotation

The main table uses plain-language classifications such as “Specific function strongly supported”, “Likely function supported by evidence”, “Protein domain detected; full function unknown”, “Conserved in phages; function unknown”, and “No reliable function identified”. Proposed function answers “what does PhageMine think it does?”, while confidence and best evidence explain how strongly and why. Machine-readable states remain in `functional_classification.tsv`. A prediction is not experimental confirmation. PhageMine does not infer lifestyle.

Batch tables include protein ID, coordinates, strand, classification, proposed function, confidence, evidence summaries, and PMF where available. Protein IDs link to detailed records rather than reducing a row to a state label alone.

## Protein retrieval

```bash
phagemine extract protein PM_000023 --run results/genome --protein-fasta
phagemine extract protein PM_000023 --run results/genome --cds-fasta
phagemine extract protein PM_000023 --run results/genome --evidence
```

Each record also links genomic context and, in discovery results, PMF membership. The internal PhageMine protein ID is the stable key connecting the report row, amino-acid FASTA, nucleotide CDS, evidence, and context.

## Important outputs

Annotation outputs include the concise `annotation.tsv`, detailed `functional_classification.tsv/json`, `evidence.json`, `candidate_ranking.tsv`, ID-only `proteins.faa`, product-labelled `annotated_proteins.faa`, `cds.fna`, `genes.gff3`, `gene_call_confidence.tsv`, `gene_calls_for_review.tsv`, `genomic_context.tsv/json`, `modules.tsv/json`, `quality_control.json`, `report.html`, `report.md`, `protein_details/`, `figures/`, `figure_data/`, and `genbank_submission/`.

`hallmark_completeness.tsv` checks whether the current evidence establishes
major capsid, portal, terminase, tail, tape-measure, replication and lysis
components. `NOT_ESTABLISHED` never means biological absence.
`annotation_review.tsv` is the short manual-curation queue combining uncertain
functions, disputed gene models and unresolved hallmark components. Detailed
classification records also state the evidence tier used for each conclusion.

Discovery outputs include `pmf_families.tsv`, `pmf_members.tsv`, `family_recurrence.tsv`, `cohort_unknown_proteome.tsv`, `conserved_neighbourhoods.tsv`, `pmfdb_validation.tsv`, `discovery_ranking.tsv`, `discovery_report.html`, `discovery_report.md`, `figures/`, `figure_data/`, and `pooled_manifest.json`.

When the INPHARED genome resource is READY, discovery also writes
`inphared_nearest_phages.tsv` as the accession-level audit trail and
`inphared_summary.tsv` as the duplicate-collapsed researcher summary.
For zero-distance Mash hits, PhageMine directly compares the query and
reference nucleotide sequences, including reverse-complement and
rotation-equivalent representations. Near-reference Mash hits still require a
formal alignment or ANI workflow; taxonomy is not inferred.

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
```

The comparison reports coordinate agreement separately from product-name
agreement. Agreement between tools is supporting computational evidence, not
experimental validation.
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
