# Methods: automated evidence-based blinded adjudication

The 689 stable units from the frozen v1.2 review sheet were assigned decisions
by a deterministic `AUTOMATED_EVIDENCE_ADJUDICATION` procedure. The procedure
read the neutralized frozen evidence packet only; it did not read the
confidential tool key or use tool identity, prior comparative performance, or
new annotation searches. Exact normalized products and the one frozen synonym
rule were accepted automatically. Compatible broader labels were classified as
`TOO_GENERAL`; ambiguous, weak, family-only, or conflicting evidence remained
`UNRESOLVABLE`. Non-informative reference products were `NOT_EVALUABLE`.

The prespecified utility weights and margins were unchanged: 1.00, 0.50, 0.25,
and 0.00 for correct, partially correct, too general, and unsupported/incorrect
decisions; unresolved observations are missing and non-evaluable observations
are excluded; margins are +0.05, -0.03, and +0.02. The neutral TSV, JSONL audit,
and openpyxl workbook were checksummed before the key was read. Unblinding then
mapped the unchanged decisions to tools. This is an automated descriptive
analysis and is not independent expert review.
