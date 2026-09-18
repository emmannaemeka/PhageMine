# FINAL CORRECTED BENCHMARK STATUS

Correction workspace: `evaluation/v1.2_tool_comparison/correction_validation_20260916/`. Frozen `final_audit_20260916` and all original evidence remain unchanged.

| Genome | PhageMine corrected | Pharokka frozen | Prokka corrected |
|---|---|---|---|
| T4 | SUCCESS | SUCCESS | SUCCESS |
| Lambda | SUCCESS | SUCCESS | SUCCESS |
| T7 | SUCCESS | SUCCESS | SUCCESS |
| T5 | SUCCESS | SUCCESS | SUCCESS |
| PhiX174 | SUCCESS | SUCCESS | SUCCESS |
| P22 | SUCCESS | SUCCESS | SUCCESS |
| Mu | SUCCESS | SUCCESS | SUCCESS |

## Root causes and resolutions

**Prokka:** The original benchmark PATH shim bypassed the installed `tbl2asn-forever` package wrapper and called the date-sensitive native binary directly. That caused exit 2 at finalization. Controlled copied-input tests show native finalization exit 1/age failure without the compatibility date and exit 0 through the supported package mechanism. Corrected runs use the installed wrapper, transparent native-exit capture, `--compliant --centre BENCH`, and a workspace-local BSD-sed adapter for the package date repair. All seven native finalization exits are 0. Remaining tbl2asn validation messages are warnings, including BadEcNumberValue and feature-quality warnings; they are preserved and reported, not suppressed. tbl2asn is required for Prokka's GenBank/Sequin submission outputs, but the GFF/FAA/FFN/TBL/TSV annotation metrics exist before finalization.

**PhiX174:** Prodigal 2.6.3 has `MIN_SINGLE_GENOME 20000`; default single mode rejects 5,386 nt. Supported meta mode runs and yields 7 secondary models. PhageMine now selects meta mode below 20,000 nt by sequence length, never by accession. PHANOTATE remains primary: 8 calls are final; secondary Prodigal models support reconciliation only and do not change final CDSs or functions. All seven corrected PhageMine runs exit 0.

## Prokka before/after explanation

See `tables/prokka_before_after.tsv`. Original pre-finalization CDS counts are recovered from logs/GFF: T4 261→261 (0.0%), Lambda 61→62 (+1.639%), T7 50→51 (+2.0%), T5 177→177 (0.0%), PhiX174 0→6 (undefined denominator), P22 61→64 (+4.918%), Mu 52→54 (+3.846%). Coordinate overlap is normalized for corrected compliant contig IDs and recorded per genome. Differences reflect Prodigal mode selection for genomes below 100,000 nt, plus corrected formatting/features; they are not interpreted as biological correctness.

## Reference-based endpoint

`tables/corrected_reference_metrics.tsv` reports exact-coordinate TP/FP/FN, precision, recall, F1, start/stop/concordance fields, named-product yield and exact product agreement on coordinate TPs. Results are independently compared with curated GenBank references; tool agreement is not ground truth. Named products are lexical output categories, not validated function. No tool is declared superior.

## Regression and provenance

The first full-suite run failed 1 test because its temporary registry was placed under the repository by the test invocation, violating the test's own default-registry invariant. The rerun used basetemp/cache outside the repository and passed 287 tests with 10 deselected. This is a test-configuration environmental failure, not a source regression. Detailed outcomes are in `validation/regression_outcomes.json`. All commands, versions, package receipt hashes, logs and output checksums are under this correction workspace. Pharokka outputs were reused from the validated frozen run. T4 PHANOTATE 294-vs-297 remains UNRESOLVED. Runtime ranking remains unsupported.

READY FOR SCIENTIFIC INTERPRETATION = YES for reference-bounded descriptive accuracy and workflow comparisons, with the stated limitations; no superiority claim is made.


Annotation warning QC found 13 native validation warning-category occurrences across corrected Prokka outputs; categories and source files are in `corrected_analysis/tables/annotation_validation_warnings.tsv`. Prokka output QC is in `corrected_analysis/validation/prokka_output_qc.tsv` (7/7 PASS). Frozen-integrity check covers 5011 files with changed count 0.


## Corrected reference endpoint (measured)

| Tool | Predicted CDS | Exact TP | FP | FN | Precision | Recall | F1 | Named products | Named-product fraction |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| PhageMine | 804 | 588 | 216 | 121 | 0.7313 | 0.8270 | 0.7762 | 439 | 0.5460 |
| Pharokka | 804 | 588 | 216 | 121 | 0.7313 | 0.8270 | 0.7762 | 398 | 0.4950 |
| Prokka | 675 | 613 | 62 | 96 | 0.9081 | 0.8622 | 0.8846 | 328 | 0.4859 |

Exact-coordinate TP requires matching accession, start, end and strand to the curated GenBank feature. Start/stop agreement and coordinate-level values are in `corrected_analysis/tables/corrected_reference_metrics.tsv`; products are compared only on coordinate TPs and exact normalized text. Prokka corrected: 675 CDS, 613 TP, precision 0.9081, recall 0.8622, F1 0.8846. PhageMine and Pharokka: 804 CDS each, 588 TP, precision 0.7313, recall 0.8270, F1 0.7762. These are reference-relative measurements on this curated panel, not universal tool rankings.

PhiX174 role separation is in `corrected_analysis/tables/phix174_caller_roles.tsv`: 8 PHANOTATE primary calls, 7 Prodigal meta-mode secondary models, 8 final CDSs; secondary models do not alter final CDSs/functions.

## Final artifact paths

- Report: `corrected_analysis/report/corrected_benchmark_report.md`
- Corrected matrix and reference metrics: `corrected_analysis/tables/corrected_completion_matrix.tsv`, `corrected_analysis/tables/corrected_reference_metrics.tsv`, `corrected_analysis/tables/corrected_aggregate_metrics.tsv`
- Prokka before/after and caller-role evidence: `corrected_analysis/tables/prokka_before_after.tsv`, `corrected_analysis/tables/phix174_caller_roles.tsv`
- Figures: `corrected_analysis/figures/corrected_completion_matrix.png`, `corrected_analysis/figures/corrected_cds.png`, `corrected_analysis/figures/corrected_named_yield.png`, `corrected_analysis/figures/corrected_reference_f1.png` (PDF and SVG companions are present)
- Validation: `corrected_analysis/validation/prokka_output_qc.tsv`, `corrected_analysis/validation/phagemine_output_qc.json`, `corrected_analysis/validation/frozen_integrity_after.json`, `validation/regression_outcomes.json`
- Artifact manifest and checksum: `artifact_manifest.json`, `artifact_manifest.sha256`
