# Database leakage and circularity controls

Before running any validation genome, determine whether the genome or nearly
identical proteins occur in PHROGs, VOGDB, Pfam, INPHARED, Swiss-Prot, or any
other source used by the engine. Record database version, release date,
checksum, sequence identity/cluster method, coverage, and the exact evidence
identifier.

Results are stratified into the following frozen operational classes:

- `EXACT_REFERENCE_OVERLAP`: exact genome sequence or protein match at ≥99%
  identity and ≥95% aligned coverage to the reference record.
- `CLOSE_HOMOLOG_OVERLAP`: no exact match, but a relevant protein match at
  90–<99% identity and ≥80% coverage, or genome ANI ≥95% with ≥80% aligned
  coverage.
- `REMOTE_HOMOLOGY`: detectable homolog/profile evidence below the close-hit
  thresholds with documented coverage and significance.
- `NO_SIGNIFICANT_REFERENCE`: no source-specific hit passes its prespecified
  significance and coverage rule; Mash distance 1.0 or zero shared hashes is
  not significance.
- `UNKNOWN_DATABASE_RELATIONSHIP`: database unavailable, incomplete, or not
  versioned sufficiently to classify.

Genome ANI and alignment coverage are recorded jointly; a single protein does
not determine genome relatedness. Results are stratified into exact/close,
remote, and no-significant-reference groups. The full panel and a
leakage-controlled subset are reported separately. Exact database presence is
provenance and is not described as independent generalization.

Nearest mathematical candidates with Mash distance 1.0, zero shared hashes,
negligible alignment, or inadequate coverage are not biologically significant
references and cannot create genus/species claims. Protein conservation,
gene-content, synteny, and module evidence are used only when their provenance
and independence are explicit.

No validation sequence is added to a database or used for threshold tuning
before the validation analysis is frozen. Any post hoc overlap analysis is
secondary and cannot change primary adjudications.
