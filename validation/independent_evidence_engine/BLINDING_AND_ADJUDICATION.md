# Blinding and adjudication

Each eligible locus receives a stable anonymous validation ID. Raw outputs,
evidence packets, reviewer decisions, and the tool key are separate artifacts.
Reviewer packets contain coordinates, reference evidence, candidate assertions,
quantitative scores, provenance, supporting and conflicting evidence, and
negative/absent-expected evidence already available from the frozen run.

The packet does not disclose tool names, tool-specific paths or filenames,
comparative summaries, expected outcomes, or post-unblinding error labels.
Presentation order is randomized with a recorded seed. A confidential mapping
is stored separately and is inaccessible to reviewers and scoring scripts until
adjudication completeness, duplicate/missing-ID checks, allowed categories,
and checksums pass.

Reviewer fields include `validation_unit_id`, `adjudication_class`,
`evidence_basis`, `reviewer_notes`, `uncertainty`, `review_required`,
`preferred_product_name`, `specificity_appropriate`, and evidence strength.
Allowed outcomes and definitions are frozen before review: `CORRECT`,
`PARTIALLY_CORRECT`, `TOO_GENERAL`, `UNSUPPORTED_SPECIFIC`, `INCORRECT`,
`UNRESOLVABLE`, and `NOT_EVALUABLE`.

Each reviewer submission is immutable. Agreement, disagreements, consensus,
and reviewer IDs are preserved independently. Unblinding is an explicit later
operation; no performance result is generated during blinded QC.

## Reviewer requirements

Primary reviewers must have documented training in bacteriophage gene/function
annotation and evidence interpretation, complete a protocol training set that
is not part of the validation panel, and pass a competency check before seeing
validation packets. They declare conflicts of interest and do not participate
in cases where they have a direct authorship or database-curation conflict.
Two independent reviewers are required for primary functional and architecture
judgments where feasible. If only one qualified reviewer is available, the
analysis is labelled single-reviewer and difficult cases receive a second
review or an `UNRESOLVABLE` status.

Disagreements remain separate until a prespecified consensus meeting. Consensus
decisions retain both original labels, the reason for resolution, reviewer
identifiers, and timestamp. Report raw agreement, Cohen's kappa for two
reviewers, or a suitable multi-rater statistic; do not replace disagreements
with a silent majority label.
