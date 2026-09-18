# v1.2 adjudication readiness report

Status: **READY FOR GENUINE BLINDED REVIEW; FINAL STATISTICS BLOCKED**.

The merged PR #14 snapshot (`4933d8f`) was audited without annotation reruns. Frozen reference and prediction checksums pass. The primary frozen counts remain: 390 evaluable reference loci, 358 PhageMine named assertions, 341 Pharokka named assertions, 292 and 38 unresolved judgments, 75 and 34 exclusive named assignments, and PhageMine per-genome yield 144/56/50/95/8/48/38.

The reviewer package contains 689 stable units (`FB0001`–`FB0689`) across exact and relaxed strata. Each unit identifies the genome, locus, coordinates, strand, reference product, neutral Prediction_A/Prediction_B products, archived raw evidence, source/provenance fields, existing scores, existing thresholds status, context/domain fields, and blank reviewer fields. The XLSX workbook has Instructions, Adjudication, and Decision_Definitions sheets with dropdown validation. The mapping is in a separate confidential file and is never included in reviewer-facing rows.

The utility definition is frozen before adjudication: CORRECT 1.00, PARTIALLY_CORRECT 0.50, TOO_GENERAL 0.25, unsupported/incorrect 0.00; UNRESOLVABLE is missing and NOT_EVALUABLE is excluded. Strict and inclusive endpoints, +0.05 functional margin, −0.03 structural margin, and +0.02 unsupported-specificity safety bound are recorded and checksummed.

The post-adjudication pipeline implements genome-level utility, strict and inclusive sensitivity endpoints, paired exact sign permutation, genome bootstrap, structural non-inferiority, unsupported-specificity safety analysis, and separate PhageMine error analysis. It refuses incomplete reviews. Synthetic complete-review tests run only in temporary directories and are labelled synthetic; no synthetic result is stored as a final result.

Tests: 30 relevant tests passed, including 7 adjudication workflow tests and the existing 23 frozen functional benchmark tests. Full final precision, recall, F1, utility, safety, agreement, and superiority conclusions remain blocked until a genuine completed review is supplied.

Next action: provide the completed blinded workbook/CSV from one or more independent reviewers, preserving case IDs, neutral evidence fields, reviewer IDs, timestamps, and decisions. Do not provide or expose the blinding key to reviewers.
