# Accuracy review: naming and gene models

This audit uses the completed fresh seven-genome searches and current-rules
decision exports. It does not run new searches, establish independent ground
truth, or demonstrate superior accuracy.

The paired PhageMine/Pharokka table contains all 804 shared predicted models
and 112 reference models missing from both tools. Across the 804 models:

| Comparison | Loci |
| --- | ---: |
| Same normalized label, including unnamed labels | 452 |
| Both unnamed, different wording | 271 |
| Different named labels | 25 |
| PhageMine unnamed, Pharokka named | 29 |
| PhageMine named, Pharokka unnamed | 27 |

The 81 differences need evidence review. Counts cover all predicted loci,
including those without an exact reference model; they are not the denominator
used in the earlier 465-named-reference-locus summary. Six of the 25 named
differences are “major capsid protein” versus “major head protein.” That is a
terminology candidate, but neither label establishes correct family-member
specificity at an individual locus. Do not resolve cases by majority vote or
reference label alone.

## Evidence and decisions

For each difference, record a stable evidence identifier and version, reviewer,
decision, and rationale in the paired table. Preserve cases without sufficient
evidence as unresolved. Review uncertainty and specificity explicitly: a domain
match does not necessarily establish whole-protein function, and a family match
does not necessarily distinguish family members. Published experimental
characterization or a traceable curated record is stronger than an annotation
copied from the same source database. Evidence from PHROGs cannot independently
validate this PHROGs-based comparison.

Review missing and changed gene models separately. Prokka's higher exact
reference-model concordance is a diagnostic signal; the reference itself still
requires assessment. Withholding a name can prevent unsupported specificity or
lose a valid annotation, so review both gains and losses before changing rules.

The scoring implementation now treats NA/N/A/NaN and DUF/UPF unknown-function
labels as non-informative. It preserves uncertainty and subunit differences.
This fixes measurement semantics, not annotation accuracy. Regression tests use
synthetic examples rather than the seven reference answers.

Rescoring the same fresh prediction tables with this policy gives:

| Tool | Named assertions at 465 named reference loci | Literal agreement | Unresolved named differences |
| --- | ---: | ---: | ---: |
| PhageMine current rules | 334 | 276 | 58 |
| Pharokka | 341 | 303 | 38 |
| Prokka | 266 | 35 | 231 |

The older summary remains an unchanged historical result under the earlier
scoring policy (335 and 268 named assertions for PhageMine and Prokka).
Predictions have not changed. Functional precision remains unset for all tools.
These counts alone cannot establish biological accuracy or rank the tools.

## Separate validation

Freeze decision rules before choosing and running the external test. Construct
a new panel with accession/version, sequence checksum, selection rationale,
reference evidence, database snapshot, and sequence-overlap assessment recorded
before inspecting tool results. Exclude the seven diagnostic genomes and assess
close relatives and database representation. A new accession alone is not proof
of an independent test. Keep development cases, rule-tuning cases, and external
test cases separate. Report reference concordance separately from independently
reviewed accuracy; do not infer unseen validation results.

Measure unsupported functional assertions, supported functional coverage,
abstention, and gene-model correctness separately. Compare tools on the same
inputs, with tool/database/search policies disclosed. Report performance per
genome and paired gains/losses; locus counts do not provide independent samples
when many loci come from the same genome. Leave functional accuracy unset while
required adjudications remain unresolved.

## Files and reproduction

- [Paired review table](results/paired-accuracy-review/paired_review.tsv)
- [Source identities and counts](results/paired-accuracy-review/manifest.json)

```bash
python evaluation/v1.3_validation/scripts/review_disagreements.py \
  --left evaluation/v1.3_validation/results/fresh-seven/inputs/PhageMine-current-rules.tsv \
  --right evaluation/v1.3_validation/results/fresh-seven/inputs/Pharokka.tsv \
  --references evaluation/v1.3_validation/results/fresh-seven/inputs/references.tsv \
  --output /path/to/new-review-directory
```

The script rejects duplicate/invalid loci and existing output directories. It
never assigns correctness automatically; the source hashes bind this diagnostic
table to the exact prediction/reference files.
