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

The primary endpoint family is **genome-level independently adjudicated
functional utility with mandatory resolution and assertion coverage**. It is a
paired estimand only when a comparator is run; it is not a superiority claim.
Exact definitions are in `STATISTICAL_ANALYSIS_PLAN.md` and are frozen before
panel outcomes are inspected.

For genome `g` and tool `t`, let `E_g` be the number of reference loci in the
primary reference tiers that are evaluable for a functional assertion. Let
`R_gt` be the subset of `E_g` with a valid resolved adjudication and let
`u_i` be the fixed utility weight for locus `i`. The primary genome-level
utility is `U_gt = sum(u_i for i in R_gt) / |R_gt|`, undefined and reported as
missing when `|R_gt|=0`. Resolution coverage is `C_gt = |R_gt| / |E_g|` and
assertion coverage is reported separately. `UNRESOLVABLE` is missing from
`U_gt`, never zero; `NOT_EVALUABLE` is outside `E_g`; missing tool output is a
failure of availability and is reported separately. The panel estimand is the
unweighted mean of genome-level `U_gt` values, with its genome-level interval,
and the paired difference when a comparator exists. No genome receives more
weight because it contains more CDSs.

The primary report always contains the pair `(mean U, mean C)` and the
unresolved fraction. A high utility with poor resolution coverage cannot be
described as broadly reliable. The protocol specifies no superiority,
non-inferiority, or numerical margin.

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

## Prespecified failure and claim controls

Formal evidence against confidence calibration is a statistically supported
reverse ordering in which HIGH has lower correctness than LOW (one-sided
ordered trend test, alpha 0.05, with direction and test frozen before review).
Failure to demonstrate a positive trend is inconclusive, not success. No
numerical unsupported-specificity safety margin is asserted without an
independent scientific basis; its estimate and interval are descriptive.

Descriptive warnings include a large unresolved fraction, loss of utility after
database-overlap control, architecture contradictions, module status based on
weak isolated evidence, and unsupported specificity concentrated in
HIGH/MODERATE calls. Review triggers include reverse confidence ordering,
repeated strong-reference contradictions, leakage-stratum divergence, or
failure of REVIEW_REQUIRED loci to have lower adjudication resolution than
ordinary loci.

Additional calibration review conditions are prespecified without arbitrary
safety margins: no meaningful improvement in independently adjudicated
reliability from LOW toward MODERATE/HIGH when category sample sizes provide
identifiable intervals; unexpectedly poor HIGH-confidence correctness;
unsupported-specificity errors concentrated in HIGH-confidence specific calls;
confidence categories whose intervals substantially overlap and provide little
discrimination; or apparently favorable calibration that disappears in the
EXACT_REFERENCE_OVERLAP/CLOSE_HOMOLOG_OVERLAP sensitivity analysis. These are
calibration-failure or major-review conditions, not automatic threshold changes
or post hoc reclassification.

Allowed conclusions are limited to observed endpoints and uncertainty:
`technically operational` requires reproducible execution; `structurally
validated`, `functionally validated`, `confidence-calibrated`, and
`biologically informative` each require their corresponding endpoint
evidence. “Improved on endpoint X” requires a paired estimate and interval for
X. “Superior,” “more accurate,” or “replaces Pharokka” is not permitted unless
a separately approved superiority analysis directly supports that exact claim.

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
