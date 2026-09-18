# Statistical analysis plan

This plan is frozen before validation data collection. It is not a power
claim and does not reuse v1.2 margins.

## Units and estimands

The primary unit is the genome. A locus is an observational unit nested within
its genome. For genome `g`, functional utility is the mean prespecified
locus utility among loci eligible for the primary reference standard; the
report includes the number evaluated, unresolved, non-evaluable, and covered.
Tool-specific coverage is named assertions divided by eligible reference loci.

The proposed utility weights, to be frozen before panel selection, are
`CORRECT=1.00`, `PARTIALLY_CORRECT=0.50`, `TOO_GENERAL=0.25`, and
`UNSUPPORTED_SPECIFIC`/`INCORRECT=0.00`. `UNRESOLVABLE` is missing and
`NOT_EVALUABLE` is excluded. These weights represent retained biological
utility with a fixed penalty for loss of specificity; they are not selected
from validation results. If protocol review changes them, the change is
versioned before any candidate is run.

Strict functional correctness counts only `CORRECT`; inclusive utility
sensitivity counts `CORRECT` plus `PARTIALLY_CORRECT` as useful. `TOO_GENERAL`
remains distinct. `UNRESOLVABLE` is missing for utility and `NOT_EVALUABLE` is
excluded. Unsupported or incorrect assertions receive zero only when the
reference is evaluable and the assertion has been adjudicated.

## Estimation

Report genome-level means, paired tool differences when a comparator is run,
and 95% confidence intervals using a paired genome bootstrap or exact
sign-flip/permutation sensitivity analysis. Locus-level intervals use cluster
resampling by genome only as a secondary descriptive analysis. No CDS is treated
as an independent biological replicate for the primary conclusion.

Binary structural and functional outcomes use paired contingency tables and
McNemar's test only when pairing and cell counts justify it. Otherwise report
exact paired results without a forced asymptotic p-value.

Confidence calibration reports, for each confidence category, observed
independent correctness, coverage, unresolved rate, and Wilson or exact 95%
intervals. Calibration is assessed descriptively with a prespecified ordered
trend statistic and a reliability plot; no threshold is changed based on the
plot. Brier score and ordinal calibration error are secondary if the reference
standard supplies a binary correctness label for the relevant loci.

Architecture, module, and novelty outcomes are summarized by genome and by
architecture stratum. Multiple secondary tests are labelled exploratory; if a
family of formal tests is declared, Holm adjustment is applied within that
family.

## Missingness and sensitivity

Missing evidence, unresolved adjudication, non-evaluable reference functions,
and unavailable databases are reported separately. Primary denominators are
never imputed. Sensitivity analyses include complete evaluable cases,
strict-correct-only scoring, inclusive useful-function scoring, and a
leakage-controlled subset. No analysis may relabel unresolved cases.

Power or sample-size calculations will be completed only after the panel
design supplies defensible estimates of genome count, clustering, and expected
reference coverage. Assumptions will be recorded rather than invented.
