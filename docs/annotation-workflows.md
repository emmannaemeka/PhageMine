# Annotation workflows in 1.3.0.dev0

These features are available on the development branch; the published 1.2 release does not include them. They improve feature coverage, interoperability and reviewability. They do not establish superior biological accuracy over another annotator.

## Annotate and assess RNA features

Install the project with `pip install .` in the tool environment described in [installation](installation.md). Run:

```bash
phagemine run genome.fasta --output results --rna-features auto
phagemine validate-run --run results
```

`auto` runs tRNAscan-SE if installed and records a skipped stage otherwise. Use `--rna-features trnascan` to require the executable, or `none` to disable it. The bacterial-mode scan exports supported simple tRNA features into RNA JSON/TSV/GFF3, the main GFF3 and the submission feature table. Ambiguous, pseudogene and intron-containing calls are recorded for specialist review. Zero predictions do not demonstrate absence. This stage is not supported on aggregate segmented runs; annotate individual segments instead.

A completed single-genome result can receive RNA features in a new directory:

```bash
phagemine rna scan --run results --output results-rna --threads 4
phagemine rna import --run results --gff rna.gff3 --source PROVIDER --version VERSION --output results-imported-rna
```

GFF3 imports accept tRNA, tmRNA, rRNA and ncRNA with matching sequence identifiers and valid coordinates. Provider name/version are mandatory. Importing features does not mean PhageMine independently detected them.

## Review additional annotation evidence

Import an existing Phold or other native GenBank annotation:

```bash
phagemine evidence-import --run results --genbank external.gbk --source Phold --version TOOL_VERSION --database-version DB_VERSION --output results-proposals
```

The nucleotide sequence must match the analyzed genome. Each proposed CDS must match an existing simple CDS in coordinates, strand and translation. Unmatched models are recorded as rejected proposals. Accepted proposals retain their source, tool/database versions and input checksum in `annotation_proposals.tsv` and `external_evidence.json`. They do not automatically replace products or change confidence. PhageMine does not run Phold or download its databases through this command.

## Apply audited review

The GUI (`phagemine gui`) includes Review / Curation and RNA Features pages. Review evidence and proposals, supply reviewer, reason and evidence reference, and save to a new result directory.

The equivalent CLI workflow is:

```bash
phagemine curate template --run results-proposals --output edits.json
# Fill reviewer and annotations in edits.json, then:
phagemine curate apply --run results-proposals --changes edits.json --output results-reviewed
```

An annotation entry looks like this (use an actual existing protein ID and supported product):

```json
{"PM_000001": {"product": "hypothetical protein", "reason": "Insufficient support for a specific product", "note": "Reviewed available evidence"}}
```

Only product, gene and note fields can be changed. Named products require an evidence reference; reviewer and reason are mandatory. Templates are tied to the current annotation-table checksum. Edits retain original evidence and automated annotation/report, update current annotation, GFF3, protein FASTA and submission files, and append a hash-chained audit trail. Coordinates and sequences remain unchanged. Existing official submission-validation outputs are invalidated when features change.

## Share and check results

```bash
phagemine validate-run --run results-reviewed
phagemine bundle create --run results-reviewed --output results-reviewed.zip
phagemine bundle verify results-reviewed.zip
phagemine profile --output runtime.json -- phagemine run genome.fasta --output timed-results
```

Validation checks artifact consistency, coordinates, translations and audit integrity. Bundle verification checks file membership and SHA-256 hashes. Neither measures annotation accuracy or NCBI acceptance. Bundles contain run artifacts and provenance, not third-party databases; retain the recorded database versions separately. Profiling reports elapsed time and sampled process-tree RSS, with caveats about shared pages and missed short peaks.

## Container and distribution

The repository includes a Dockerfile and a container smoke workflow:

```bash
docker build -t phagemine:dev .
docker run --rm -v "$PWD:/data" phagemine:dev run /data/genome.fasta --output /data/results --rna-features trnascan
```

The base-image tag and several external dependencies are not fully locked. Resource databases are configured separately. The local Bioconda recipe under `packaging/bioconda/recipe` is a development recipe; it has not been submitted to or published by Bioconda. Consult CI results before treating a container build as verified.

## Database health and reusable searches

```bash
phagemine database-health --checksums --output database-health.json
phagemine database-health --expected-versions expected-versions.json --max-age-days 180
phagemine run genome.fasta --output first-run --evidence-cache /path/to/cache
phagemine run genome.fasta --output updated-run --evidence-cache /path/to/cache
```

Expected versions are an explicit JSON map such as `{"PHROGS": "YOUR_INSTALLED_RELEASE"}`. Checks run offline: they detect missing resources, broken preparation/index files, version mismatches and ambiguous installations. Recorded release or installation dates support maintenance warnings; registration time and file modification time are not treated as database freshness. Missing optional resources allow core analysis with a limited-evidence status. Version expectations with no matching registration fail the health check. Every single-genome run saves `database_health.json` before prediction.

Persistent reuse is opt-in for single-genome evidence searches. Only successful real searches and successful zero-hit results are stored. Sequence, database/index/annotation bytes, executable bytes, versions, adapter parameters and package implementation must match. Corrupt entries are recomputed; changed inputs trigger fresh searches. Reading and hashing large databases takes time. Gene prediction and RNA stages still run; manual edits use the existing curation command. Cache provenance labels old search commands as historical. Keep the cache separate from result directories when exporting bundles.

## Evidence conflicts

`annotation_conflicts.tsv` and `.json` expose conflicts already flagged by the fusion engine, plus differences between a current named product and an imported external proposal. A difference in wording does not prove incompatible function. The GUI can filter these records and display the underlying evidence during review. Saving an edit refreshes the conflict report without automatically raising confidence or resolving remaining alternative proposals.

## Curated-reference benchmark

```bash
phagemine benchmark-curated --predictions predictions.tsv --references references.tsv --output benchmark-results
phagemine benchmark-curated --predictions predictions.tsv --references references.tsv --adjudications reviewed-decisions.tsv --output benchmark-reviewed
```

Both tables require `accession`, `start`, `end`, `strand` and `product`. References also require a nonempty `curated_by` and `reference_evidence` on every row. Accessions and coordinates must refer to the same representation; matching is exact and does not silently pool different genomes. Missing exact models count as abstentions at informative reference loci, and unmatched predicted models are reported separately.

Reports include model precision/recall, product precision, coverage, abstention counts and Wilson intervals. Exact normalized labels are scored automatically. Other named-label differences enter `review_queue.tsv`; precision/recall scores affected by unresolved decisions remain unset. Reviewers can fill category, reviewer and rationale and pass the completed queue as adjudications. Optional `--synonyms` accepts a TSV with `term_a`, `term_b`, `rule_id` and `rationale`. Decisions and source tables are checksummed. Provide independently curated references; the software cannot establish reference independence or quality from a curator name. This workflow does not itself establish superiority over another tool.
