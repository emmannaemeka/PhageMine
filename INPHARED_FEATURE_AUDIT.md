# INPHARED feature audit — PhageMine v1.2

Audit scope: source, CLI, resource manager/installer, pipeline, discovery
workflow, tests, documentation, examples, and changelog. The seven-genome
structural benchmark was not rerun.

## Classification

**A — Fully implemented and currently callable.** INPHARED is integrated into
the single-genome `annotate` and `run` workflows and into batch `discover` and
`both` workflows. It is automatically consumed only after operational resource
validation; otherwise the workflow records a safe `SKIPPED` result.

## Implementation

- `src/phagemine/inphared.py` implements resource preparation, Mash screening,
  exact sequence confirmation for zero-distance matches, and optional
  bidirectional BLASTN length-normalized similarity.
- `src/phagemine/resources.py` validates canonical INPHARED files, checksums,
  manifests, and the required `mash` executable. It reconstructs runtime paths
  from the validated resource root rather than stale provenance paths.
- `src/phagemine/database_installer.py` downloads the pinned 7 April 2026
  INPHARED genome and metadata artifacts, prepares a per-genome Mash sketch,
  writes QC and manifests, and registers the resource only after validation.
  The same release supplies PMFDB's protein reference conversion, which is a
  separate protein-family comparison.
- `src/phagemine/pipeline.py` invokes the validated comparison for `annotate`
  and `run`; `src/phagemine/batch.py` resolves validated comparative resources
  for batch discovery/both; `src/phagemine/discovery.py` writes comparative
  outputs and inserts an idempotent report section.

## Invocation and setup

There is no separate `inphared` analysis subcommand. Register or install the
resource, then run the normal workflow:

```bash
phagemine databases install inphared
phagemine doctor
phagemine run genome.fasta --output results/genome
phagemine batch genomes/ --output results/cohort --mode discover --evidence core
```

The installer downloads the pinned INPHARED release `2026-04-07` from the
Millard Lab S3 distribution, verifies published MD5 values, converts the
multi-FASTA and metadata table, and runs `mash sketch -i`. The prepared root
must contain `reference_phage_genomes.fna`, `genome_metadata.tsv`,
`genome_qc.tsv`, `genome_manifest.json`, and `inphared.msh`; `mash` must be
operationally available.

## Algorithms and scientific meaning

Mash `dist` screens each query against the per-genome INPHARED sketch and ranks
the nearest references by Mash distance, p-value, and matching hashes. Mash is
screening evidence only; PhageMine does not convert distance into ANI.

For zero-distance hits, the implementation can confirm identical, reverse-
complement, or rotation-equivalent nucleotide sequence when the reference
FASTA is available. It can then run BLASTN in both directions and report
PhageMine's length-normalized intergenomic similarity and query/reference
coverage. This is not VIRIDIC or ANI.

Metadata returned from INPHARED includes reference accession, description,
host genus, phage genus, subfamily, and family. These are reference metadata,
not an assignment of the query. A 95% species working boundary is labelled as
requiring ICTV review; no generic genus boundary is configured, and partial or
poorly aligned queries are explicitly ineligible for numerical taxonomy. The
code does not make formal ICTV taxonomic assignments.

## Outputs

The comparison directory contains:

- `inphared_nearest_phages.tsv` and `.json`: accession-level ranked hits;
- `inphared_summary.tsv`: duplicate-collapsed nearest-reference summaries;
- run manifest `comparative_analysis.inphared` and report section.

Columns include `sample_id`, `rank`, `reference_accession`, `mash_distance`,
`p_value`, `matching_hashes`, `reference_description`, `host_genus`,
`phage_genus`, `phage_subfamily`, `phage_family`,
`intergenomic_similarity_percent`, query/reference aligned percentages,
`genome_length_ratio`, `similarity_method`, species/genus thresholds and
source, `taxonomy_eligible`, `comparison_scope`,
`taxonomic_interpretation`, sequence confirmation, and interpretation.

When the resource is unavailable, incomplete, checksum-invalid, or missing
Mash, PhageMine records `status=SKIPPED` and the reason. It does not create
zero-hit biological claims. Mash or BLASTN execution failures are propagated as
workflow errors rather than silently converted to success.

## Tests and release implications

The existing INPHARED/resource tests pass (`25 passed` for
`tests/test_inphared.py` and `tests/test_database_installer.py`); with the
small-genome regression and release-metadata checks, the validated targeted
subset was `32 passed`. They cover ranked output, duplicate accession collapse,
sequence confirmation, BLASTN fallback, ICTV-boundary abstention, canonical
path resolution, checksum/resource failures, installer preparation, and report
idempotency.

The completed seven-phage structural benchmark did not validate INPHARED
nearest-neighbour or taxonomic performance. Public documentation must present
INPHARED as comparative genomic context and reference metadata, with cautious
similarity interpretation, rather than as a validated taxonomic classifier.

## Release blockers

No algorithmic blocker was found. Before release, retain the pinned source
release and checksums in user run manifests, keep large INPHARED files out of
GitHub, and ensure the release environment supplies a working `mash` and
optionally `blastn` for intergenomic similarity.
