# Database leakage and circularity controls

Before running any validation genome, determine whether the genome or nearly
identical proteins occur in PHROGs, VOGDB, Pfam, INPHARED, Swiss-Prot, or any
other source used by the engine. Record database version, release date,
checksum, sequence identity/cluster method, coverage, and the exact evidence
identifier.

Results are stratified into exact/reference-sequence present, close cluster,
remote homolog, and no significant reference. The full panel and a
leakage-controlled subset are reported separately. Exact database presence is
reported as provenance and is not described as independent generalization.

Nearest mathematical candidates with Mash distance 1.0, zero shared hashes,
negligible alignment, or inadequate coverage are not biologically significant
references and cannot create genus/species claims. Protein conservation,
gene-content, synteny, and module evidence are used only when their provenance
and independence are explicit.

No validation sequence is added to a database or used for threshold tuning
before the validation analysis is frozen. Any post hoc overlap analysis is
secondary and cannot change primary adjudications.
