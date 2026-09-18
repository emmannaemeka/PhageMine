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
