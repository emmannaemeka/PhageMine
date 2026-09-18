# PhageMine v1.2 benchmark final report

## Status

This report closes the frozen v1.2 benchmark and its separate
`AUTOMATED_EVIDENCE_ADJUDICATION` analysis track. The original human-review
workflow remains preserved and blank for future independent validation. No
annotation was rerun and no frozen prediction, reference, coordinate, matching,
synonym, scoring, or utility input was changed.

The 689 stable blinded units were processed using only their frozen neutral
evidence packets. Automated decisions were frozen and checksummed before the
confidential tool key was accessed.

## Automated adjudication

| Outcome | Count |
|---|---:|
| CORRECT | 0 |
| PARTIALLY_CORRECT | 161 |
| TOO_GENERAL | 52 |
| UNSUPPORTED/INCORRECT | 0 |
| UNRESOLVABLE | 476 |
| NOT_EVALUABLE | 0 |
| **Total** | **689** |

Resolved/scored decisions were 213/689 (30.9144%); unresolved decisions were
476/689 (69.0856%). Confidence counts were HIGH 38, MODERATE 175, and LOW 476.

## Utility

The frozen utility weights remain CORRECT=1.00, PARTIALLY_CORRECT=0.50,
TOO_GENERAL=0.25, and UNSUPPORTED/INCORRECT=0.00. UNRESOLVABLE observations are
missing and NOT_EVALUABLE observations are excluded.

The observed missing-aware bounds were:

| Tool | Observed lower bound | Observed upper bound |
|---|---:|---:|
| PhageMine | 0.107692 | 0.869231 |
| Pharokka | 0.004487 | 0.986538 |

The observed lower-bound difference (PhageMine minus Pharokka) was **+0.103205**.
It numerically exceeds the prespecified +0.05 margin, but the high unresolved
fraction prevents a valid confidence-interval-based superiority inference. No
confidence interval was manufactured, and unresolved cases were neither
imputed nor treated as zero.

The automated adjudicator assigned no cases to UNSUPPORTED/INCORRECT. This is a
descriptive property of this automated rule set, not evidence that either tool
has zero true unsupported-specificity errors.

## Structural result

The frozen structural benchmark is retained unchanged. PhageMine and Pharokka
have the published identical structural results; this final report does not
rerun structural annotation or infer a new statistical interval.

## Interpretation

The result is an automated evidence-based blinded adjudication, not independent
human expert adjudication. It should be used as a reproducible descriptive
analysis and readiness record for later independent human validation, not as a
claim of comparative superiority.
