# Independent evidence-engine validation protocol

## Scope and separation

This prospective study begins from protected-main commit
`318e963db1288c400925e3eaaf02e512f774949a`. It does not modify the engine,
thresholds, PHANOTATE policy, Pyrodigal-gv reconciliation, v1.2 benchmark,
or v1.2 adjudication. The development cases `k141_1259522`, `k141_247505`,
`k141_371497`, `k141_731777`, and `k141_516427`, and the seven prior reference
phages T4, Lambda, T7, T5, PhiX174, P22, and Mu, are permanently excluded.

## Questions

1. Does the PHANOTATE-final, Pyrodigal-gv-observational, biological-support
   framework distinguish exact, boundary-tolerant, short, overlapping, and
   conflicting CDS evidence without treating caller agreement as truth?
2. Do functional classes separate specific correctness, useful generality,
   domain-only/family evidence, unsupported specificity, and abstention?
3. Are HIGH/MODERATE/LOW/UNRESOLVED structural and functional confidence
   categories ordered by independently verified reliability?
4. Do module aggregates preserve locus provenance and avoid promoting weak
   repeated hits into module certainty?
5. Are architecture hypotheses and architecture-specific hallmarks coherent,
   including appropriate non-applicability for filamentous phages?
6. Does comparative analysis distinguish nearest candidates from significant
   references and avoid automatic novelty/taxonomy claims?
7. Do LOW, UNRESOLVED, CONFLICTING, and REVIEW_REQUIRED cases identify work
   that is genuinely difficult for an independent reviewer?

## Primary endpoint

The primary endpoint is **genome-level independently adjudicated functional
utility among eligible, evaluable reference loci, reported together with
coverage and unresolved fraction**. Utility is calculated as the mean of
prespecified adjudication weights for each genome; unresolved and non-evaluable
observations remain missing/excluded under the reference-standard rules and
are never silently scored as incorrect. The primary report contains one
estimate per genome and a confidence interval that respects genome clustering.
The protocol does not specify superiority, non-inferiority, or a margin.

This endpoint is selected because it rewards defensible information while
making abstention and denominator changes visible. Exact definitions and
sensitivity estimands are frozen in `STATISTICAL_ANALYSIS_PLAN.md` before
panel outcomes are inspected.

## Secondary endpoints

- CDS exact-coordinate, boundary-tolerant, start/stop, short-CDS, and
  overlap-aware agreement with the independent reference.
- Functional correctness, specificity, unsupported-specificity rate,
  excessive-generalization rate, unresolved rate, and evidence-class recovery.
- Calibration and discrimination of confidence categories.
- Module-level status accuracy and provenance completeness.
- Architecture classification and architecture-specific hallmark agreement.
- Comparative significance/conservatism, including false novelty claims.
- Reviewer time, review-required fraction, and inter-reviewer agreement.

## Reference and adjudication

Reference records are assembled under the hierarchy in `REFERENCE_STANDARD.md`.
Product strings alone are insufficient. Reviewers receive stable blinded unit
IDs and evidence packets containing supporting, conflicting, weak, rejected,
and absent-expected evidence. Tool identity is held in a separate key until
adjudication and QC are frozen.

## Analysis guardrails

Loci are nested in genomes; genome summaries are the primary inferential unit.
No validation genome may be used for threshold tuning. Any exploratory
analysis, database-overlap subset, or comparator result is labelled secondary
and descriptive unless its estimand is explicitly locked in advance.
