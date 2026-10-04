# Annotation benchmark and accuracy review

The seven-genome benchmark and review of all 81 PhageMine/Pharokka differences
are complete. The results do **not** establish that PhageMine is more accurate
overall. Gene-model concordance, supported functional coverage and unsupported
specificity must be assessed separately.

## Fresh benchmark

[All seven jobs passed](https://github.com/emmannaemeka/PhageMine/actions/runs/37182016577).
PhageMine, Pharokka 1.10.1 and Prokka 1.14.6 generated new predictions from
identical nucleotide-only inputs. PhageMine and Pharokka used PHROGs v4 from
the same database distribution; PhageMine used PHROGs only. Prokka used its
packaged viral databases. This is a diagnostic comparison, not a held-out
accuracy trial.

| Tool | Predicted CDS | Exact reference models / 700 | Named functions at 465 named reference loci | Literal name agreement |
| --- | ---: | ---: | ---: | ---: |
| PhageMine | 804 | 588 | 334 | 276 |
| Pharokka | 804 | 588 | 341 | 303 |
| Prokka | 678 | 616 | 266 | 35 |

The current scoring policy excludes missing descriptions and unknown-function
family labels from named functions. Literal agreement measures wording, not
biological accuracy. Prokka's terminology differs substantially. Compound CDS
are excluded. Reference annotations and database overlap are not independent
ground truth.

Fusion rules 1.11 were applied to the fresh evidence. All seven exports passed
integrity checks; raw evidence and gene coordinates stayed unchanged. There
were **zero product changes** relative to this fresh rules-1.10 run because the
new safeguards concern reviewed whole-protein/domain combinations absent from
the PHROGs-only configuration. This does not test the full Swiss-Prot/Pfam arm.

## Completed evidence review

All 81 differences were assessed. **58** predicted sequences exactly match
versioned reviewed UniProt records; **23** have no exact reviewed match in the
snapshot and remain unresolved. No decision was inferred from a reference name
or another tool's vote alone.

| Evidence-review outcome | PhageMine | Pharokka |
| --- | ---: | ---: |
| Named role consistent with the reviewed record | 20 | 18 |
| Compatible broad role | 4 | 7 |
| Function withheld | 20 | 21 |
| Specificity not established by the record | 6 | 5 |
| Assertion not established by the record | 7 | 4 |
| Function conflicts with the reviewed record | 1 | 3 |
| No exact reviewed match | 23 | 23 |

These are selected disagreements, not an unbiased sample for accuracy
percentages or tool rankings. Reviewed records may themselves contain inferred
or uncertain functions. An uncharacterized record does not prove a predicted
function false. The review was performed by Codex, not an independent expert.
Strict functional accuracy remains unset rather than scoring unresolved cases
as correct or incorrect.

The review identifies missed functions in PhageMine as well as useful functions
that Pharokka withheld. Naming-only differences, such as head/capsid maturation
protease, must be distinguished from genuine family-member confusion. Broad
head-morphogenesis labels can also conceal a more useful assembly or scaffolding
role. These findings support adding discriminating curated evidence, rather
than optimizing label counts against the reference.

## Retained correction and naming safeguards

The T4 protein at 107323–108606 on the plus strand matches reviewed **P19896**
exactly. Its supported product is **capsid vertex protein**, which forms pentons,
distinct from the major capsid protein forming hexamers. The user-approved edit
is retained with a source-bound curation file and audit trail. Its regenerated
outputs pass integrity checks. The original automatic benchmark is unchanged.

The naming engine now preserves a supported reviewed whole-protein name when a
generic domain rule is also present, including uncertainty words. Contradictory
reviewed names remain unresolved. Synthetic regression tests cover these rules;
the engine contains no T4 locus-specific naming override.

## Additional panel

Three dairy Lactococcus phages—c2, SK1 and Tuc2009—were selected before viewing
their tool results. Their reference versions, sequence hashes and selection
record are frozen in [external_panel](external_panel/). They were not used to
adjust the naming rules. Database overlap has not been excluded, so this panel
tests transfer outside the original genomes, not independent functional truth.
The expanded fresh comparison workflow generates new predictions for them.

## Review files

- [Benchmark summary and manifest](results/accuracy-review/benchmark_manifest.json)
- [All 81 evidence reviews](results/accuracy-review/reviewed_cases.tsv)
- [Review counts and interpretation](results/accuracy-review/review_summary.json)
- [Versioned reviewed records](results/accuracy-review/reviewed_records.json)
- [Capsid vertex correction](results/accuracy-review/capsid_vertex_curation.json)
- [Retained corrected result](results/accuracy-review/capsid_vertex_result.json)

The review validator checks prediction wording, exact protein sequences,
reviewed-record hashes and complete case coverage. Run
`scripts/complete_accuracy_review.py --help` for reproduction arguments. Fresh
workflow artifact identifiers and ZIP hashes are retained with the data.
