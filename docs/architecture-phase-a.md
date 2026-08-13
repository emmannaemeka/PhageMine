# Phase A architecture review

## Preserved MVP components

The MVP already provides structured evidence records, independent functional-confidence/biological-interest/evidence-diversity outputs, discovery ranking and reports, QC, provenance, GenBank pre-submission packaging, metadata/readiness states, optional `table2asn`, and tests.

## This phase

The direct `simple_orf_demo` call is replaced in production flow by a `GenePredictor` interface and `PHANOTATEPredictor`. PHANOTATE remains an external, optional executable: it is not installed or bundled. The demonstration predictor remains available only through explicit selection for fixtures and tests. Predictor name, version, executable, and parameters enter the run manifest.

## Deferred phases

Lifecycle prediction, curated core signatures, evidence packs/database manager, real similarity/HMM integrations, GUI, `table2asn` real-world validation, and benchmarking remain deferred. No biological database is downloaded or committed.
