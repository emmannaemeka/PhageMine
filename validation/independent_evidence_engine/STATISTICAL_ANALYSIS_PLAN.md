# Statistical analysis plan

This plan is frozen before validation data collection. It is not a power
claim and does not reuse v1.2 margins.

## Units and estimands

The primary unit is the genome. A locus is an observational unit nested within
its genome. For genome `g`, functional utility is the mean prespecified
locus utility among loci eligible for the primary reference standard; the
report includes the number evaluated, unresolved, non-evaluable, and covered.
Tool-specific coverage is named assertions divided by eligible reference loci.
The primary reference set is TIER_1_EXPERIMENTALLY_ESTABLISHED,
TIER_2_EXPERT_CURATED, and TIER_3_STRONGLY_ORTHOLOGY_SUPPORTED; TIER_4 is
reserved for the expanded-tier sensitivity analysis.

The proposed utility weights, to be frozen before panel selection, are
`CORRECT=1.00`, `PARTIALLY_CORRECT=0.50`, `TOO_GENERAL=0.25`, and
`UNSUPPORTED_SPECIFIC`/`INCORRECT=0.00`. `UNRESOLVABLE` is missing and
`NOT_EVALUABLE` is excluded. The 0.50 value treats a materially useful but
incomplete assertion as half of a fully supported assertion; 0.25 represents
compatibility with substantial loss of informative specificity; and 0.00
prevents unsupported or conflicting specificity from earning utility. These
are ordinal-to-numeric assumptions, not natural biological distances. They are
therefore a prespecified composite and must be accompanied by strict,
inclusive, ordinal-distribution, and coverage-versus-correctness sensitivities.
They are not selected from validation results. If protocol review changes them,
the change is versioned before any candidate is run.

Strict functional correctness counts only `CORRECT`; inclusive utility
sensitivity counts `CORRECT` plus `PARTIALLY_CORRECT` as useful. `TOO_GENERAL`
remains distinct. `UNRESOLVABLE` is missing for utility and `NOT_EVALUABLE` is
excluded. Unsupported or incorrect assertions receive zero only when the
reference is evaluable and the assertion has been adjudicated.

## Estimation

Report one `U_gt` and `C_gt` per genome, then unweighted panel means, paired
tool differences when a comparator is run, and 95% confidence intervals using
a paired genome bootstrap or exact sign-flip/permutation sensitivity analysis.
Locus-level intervals use cluster resampling by genome only as a secondary
descriptive analysis. No CDS is treated as an independent biological replicate
for the primary conclusion. A genome with no resolved eligible loci remains
missing for utility but contributes zero resolved coverage and is retained in
the availability report.

Binary structural and functional outcomes use paired contingency tables and
McNemar's test only when pairing and cell counts justify it. Otherwise report
exact paired results without a forced asymptotic p-value.

Confidence calibration reports, for each confidence category, observed
independent correctness, coverage, unresolved rate, and Wilson or exact 95%
intervals. Calibration is assessed descriptively with a prespecified ordered
trend statistic and a reliability plot; no threshold is changed based on the
plot. Brier score and ordinal calibration error are secondary if the reference
standard supplies a binary correctness label for the relevant loci.
HIGH, MODERATE, and LOW are tested for monotonic reliability; UNRESOLVED is
reported as a separate abstention/triage category and is never placed on that
confidence scale.

Calibration review also examines whether reliability improves meaningfully from
LOW toward MODERATE/HIGH, whether HIGH-confidence correctness is poor, whether
unsupported specificity is concentrated in HIGH calls, whether category
intervals provide useful discrimination, and whether apparent calibration
vanishes after exact/close database-overlap control. These conditions use the
reported estimates, confidence intervals, ordered analyses, and documented
biological review; no arbitrary numerical safety margin is introduced.

Architecture, module, and novelty outcomes are summarized by genome and by
architecture stratum. Multiple secondary tests are labelled exploratory; if a
family of formal tests is declared, Holm adjustment is applied within that
family.

## Missingness and sensitivity

Missing evidence, unresolved adjudication, non-evaluable reference functions,
and unavailable databases are reported separately. Primary denominators are
never imputed. Sensitivity analyses include CORRECT-only strict correctness;
CORRECT plus PARTIALLY_CORRECT useful-function correctness; ordinal
adjudication distributions; resolved-annotation correctness;
coverage-versus-correctness plots; primary versus expanded reference tiers;
database-overlap strata; architecture and genome-size strata;
clustering-threshold sensitivity; exclusion of exact-reference-overlap cases;
unresolved as a separate outcome; and resolved-only analyses with mandatory
coverage. No analysis may relabel unresolved cases or choose a favorable
sensitivity definition after outcomes.

Power or sample-size calculations will be completed only after the panel
design supplies defensible estimates of genome count, clustering, and expected
reference coverage. Assumptions will be recorded rather than invented.

## Primary and secondary hypothesis families

The primary analysis estimates genome-level utility and resolution coverage,
not superiority. Key secondary hypotheses are (1) confidence correctness is
non-decreasing from LOW to MODERATE to HIGH, (2) unsupported specificity is
described separately from utility, and (3) architecture/module statuses agree
with independent evidence. These key secondary tests use Holm adjustment
within the declared family when formal p-values are reported. Structural,
overlap, taxonomy, reviewer-burden, and detailed error analyses are
exploratory unless explicitly promoted before data collection.
