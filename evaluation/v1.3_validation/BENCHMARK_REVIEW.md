# Annotation benchmark and evidence review

The [fresh full-system benchmark](https://github.com/emmannaemeka/PhageMine/actions/runs/37186254566)
completed on all ten genomes with **all six resources**. The comparisons below
use PhageMine's final adjudicated products, not individual PHROGs hits.
PHROGs, Swiss-Prot, Pfam and VOGDB provide annotation evidence; PMFDB and
INPHARED provide comparative context rather than functional-name votes.
Every required adapter and both comparisons completed. The downloaded artifact
and all 5,202 internal file checksums were verified. All ten nucleotide inputs
are identical to the earlier Pharokka and Prokka comparison inputs.

## Full-system findings

| Configuration | Predicted CDS | Exact reference models / 851 | Named assertions at 537 informative reference loci | Exact or equivalent names | Unresolved named differences |
| --- | ---: | ---: | ---: | ---: | ---: |
| PhageMine, all resources, fresh rules 1.12 | 978 | 719 | 395 | 114 | 281 |
| PhageMine, all resources, corrected rules 1.13 | 978 | 719 | 409 | 137 | 272 |
| Pharokka 1.10.1, earlier matched-input run | 977 | 719 | 407 | 362 | 45 |
| Prokka 1.14.6, earlier matched-input run | 824 | 747 | 273 | 39 | 234 |

These measure reference concordance, **not independently established functional
accuracy**. Reviewed Swiss-Prot descriptions often use different terminology
from RefSeq and PHROGs. A change from “exonuclease” to “exodeoxyribonuclease”, for
example, requires semantic review rather than a wrong-name verdict. The scorer
uses the existing frozen equivalences; no new synonyms were added to improve
this result. Reference names can also be wrong: the T4 vertex protein illustrates
why maximising wording agreement is not an appropriate objective.

The full run automatically selects **capsid vertex protein** for T4
107323–108606 (+). Its reviewed Swiss-Prot P19896 alignment has 100% identity
and complete query/reference coverage. The selected name is distinct from the
PHROGs/VOGDB major-capsid label, and requires no manual override.

Two general decision-rule defects were exposed and corrected:

- Extended unknown descriptions, such as “uncharacterized 7.3 kDa protein in
  an intergenic region”, could previously be treated as functional names and
  displace informative candidates. Reviewed status does not establish function.
- Descriptions consisting only of `protein` and the database record's gene
  symbol, such as “Protein rIIA”, could replace an informative role. The record
  remains evidence, but that description supplies no transferable function.

Rules 1.13 changed 184 product/gene rows when applied to the fresh all-database
evidence. All ten regenerated exports passed integrity checks; evidence objects,
gene coordinates and sequences were retained. This is explicitly a decision-rule
reanalysis of the newly generated full-system evidence, not another database
search. Informative alternatives remain computational hypotheses, not functions
proved by discarding an unknown label. No manual product edits enter the scores.
There are 312 reviewed-anchor selections, and 24 proteins retain conflicting
evidence. Confidence remains rule based and uncalibrated. Local validation:
417 tests passed, two skipped, 12 external/integration tests deselected.

## Why Prokka matches more gene boundaries

Prokka used Prodigal 2.6.3 (`-c -m -g 11 -p single`), whereas PhageMine used
PHANOTATE. The full database run retained all 978 PhageMine gene models; adding
annotation databases therefore did not address gene-boundary differences.

Both tools match 684 reference models. Prokka alone matches 63, and PhageMine
alone matches 35. Of Prokka's 63, **53 share the PhageMine stop position but have
a different start**. Of PhageMine's 35, 13 share a Prokka stop position. This
identifies start-site selection as the main source of Prokka's net advantage
on this panel. It does not prove every reference start correct or justify
replacing all PHANOTATE predictions. The next gene-model experiment should
compare alternative caller starts under the same annotation evidence and
review which start is supported, retaining small and overlapping phage genes.

The full annotation commands took a median 229.7 seconds on a two-thread GitHub
runner, excluding database installation. This includes additional searches and
comparative work absent from the earlier commands; it is not a controlled speed
comparison. Database versions are retained: Swiss-Prot 2026_03, Pfam 38.2,
VOGDB 235, PHROGs v4/Pharokka 1.11.0 and PMFDB/INPHARED 2026-04-07. The earlier
PHROGs-only comparison used the Pharokka 1.8.0 distribution, so changes cannot
be attributed exclusively to adding databases.

## Earlier PHROGs-only diagnostic baseline

[All ten fresh comparisons succeeded](https://github.com/emmannaemeka/PhageMine/actions/runs/37183614867).
PhageMine, Pharokka 1.10.1 and Prokka 1.14.6 received identical nucleotide-only
inputs. Reference gene names and product descriptions were used only for
scoring. PhageMine and Pharokka shared PHROGs v4; Prokka used its packaged
viral databases. Eleven compound reference CDS were excluded.

| Tool | Predicted CDS | Exact reference models / 851 | Named assertions at 537 informative reference loci | Exact or equivalent names |
| --- | ---: | ---: | ---: | ---: |
| PhageMine, PHROGs only | 978 | 719 | 394 | 346 |
| Pharokka | 977 | 719 | 407 | 362 |
| Prokka | 824 | 747 | 273 | 39 |

These are reference-concordance counts, not functional accuracy percentages.
Prokka uses substantially different terminology; low wording agreement does
not prove its functions wrong. Matching gene boundaries does not independently
validate the biological gene model. The scoring policy excludes missing
labels, unknown-function families and identifier-only descriptions such as
“Orf80”; informative role descriptions containing identifiers remain eligible.
Whole-term equivalence rules preserve uncertainty, subunits and specificity.
“Capsid vertex protein” is never equated with “major capsid protein”.

The original seven genomes contribute 700 simple reference CDS and 451
informative reference loci. Their exact-or-equivalent name counts are 290,
303 and 38 for PhageMine, Pharokka and Prokka respectively.

Three dairy Lactococcus phages—c2, SK1 and Tuc2009—were selected before
inspecting their outputs and were not used to tune naming rules.

| Tool | Predicted CDS | Exact reference models / 151 | Named assertions at 86 informative reference loci | Exact or equivalent names |
| --- | ---: | ---: | ---: | ---: |
| PhageMine, PHROGs only | 174 | 131 | 60 | 56 |
| Pharokka | 173 | 131 | 66 | 59 |
| Prokka | 146 | 131 | 14 | 1 |

PhageMine and Pharokka share the original-panel coordinates. On the external
panel, PhageMine adds one hypothetical CDS on SK1 at 27591–27689 on the minus
strand. All ten fresh PhageMine exports were reclassified with fusion rules
1.12: raw evidence and coordinates were preserved, no automatic product
changed, and every export passed artifact-consistency checks. The original raw
searches used rules 1.11. New naming safeguards therefore have no measured
benefit in this PHROGs-only configuration.

## Earlier selected-disagreement evidence review

All 81 selected PHROGs-only PhageMine/Pharokka differences from the original
panel were reviewed. These verdicts do not score the new full-system predictions. Fifty-eight predicted sequences exactly match versioned reviewed
UniProt records; 23 have no exact reviewed match in the snapshot and remain
unresolved. This is a Codex review of selected disagreements, not independent
expert adjudication or an unbiased precision sample.

| Review outcome | PhageMine | Pharokka |
| --- | ---: | ---: |
| Named role consistent with reviewed record | 20 | 18 |
| Compatible broad role | 4 | 7 |
| Function withheld | 20 | 21 |
| Specificity not established | 6 | 5 |
| Assertion not established | 7 | 4 |
| Function conflicts with record | 1 | 3 |
| No exact reviewed match | 23 | 23 |

Reviewed records may contain inferred or uncertain functions. An uncharacterized
record does not prove a proposed function false. Functional accuracy remains
unset; unresolved cases are not silently counted as correct or incorrect.

## Informative names and the retained correction

The T4 protein at 107323–108606 on the plus strand exactly matches reviewed
P19896. Its supported product is **capsid vertex protein**, which forms pentons,
distinct from the major capsid protein forming hexamers. The user-approved
correction from the earlier baseline is retained with a source-bound curation
file and audit trail. The full-system automatic result is separate. Its
exports pass integrity checks. This manual edit is excluded from automatic
benchmark scores.

The naming rules retain informative supported roles and uncertainty words.
Identifier-only labels and DUF/UPF families do not establish functions. A domain
label cannot replace a supported reviewed whole-protein description; contradictory
reviewed names remain unresolved. There is no T4 locus-specific naming override.

## What remains unproven

The full-system run establishes that all evidence paths executed and identifies
concrete naming defects and an automatically supported vertex-protein correction.
It does not establish overall functional-accuracy superiority. Coarse family assignments can confuse distinct
structural proteins; conservative partial-match gates can withhold supported
functions; some selected descriptions can assert more specificity than their
evidence establishes. Gene-model decisions and confidence rules also need
independent calibration. More named proteins alone do not demonstrate improvement.

The panel is small, reference annotations are imperfect and database overlap
has not been excluded. The external panel measures transfer outside the original
genomes, not independent functional truth. The much larger full-system semantic review remains unfinished. The lower
wording agreement must not be described as a measured biological accuracy loss
or hidden by broad synonym substitutions.

## Data

- [Full-system comparison counts](results/full-system/panel_summary.tsv)
- [Full-system review provenance and limitations](results/full-system/review_manifest.json)
- [All adjudicated names, rules 1.13](results/full-system/all_product_names.tsv)
- [Rule changes from the fresh full output](results/full-system/decision_changes.tsv)
- [Gene-model differences](results/full-system/gene_model_differences.tsv)
- [All-resource execution checks](results/full-system/source_audits.json)
- [Preserved-evidence integrity audit](results/full-system/reclassification_audit.json)
- [Full-system checksums](results/full-system/SHA256SUMS)

Earlier diagnostic data:

- [Panel counts and denominators](results/accuracy-review/panel_summary.tsv)
- [Benchmark provenance](results/accuracy-review/benchmark_manifest.json)
- [All 978 names, including the separately marked correction](results/accuracy-review/all_product_names.tsv)
- [All 81 evidence reviews](results/accuracy-review/reviewed_cases.tsv)
- [Review interpretation](results/accuracy-review/review_summary.json)
- [Versioned reviewed records](results/accuracy-review/reviewed_records.json)
- [Capsid vertex correction](results/accuracy-review/capsid_vertex_curation.json)
- [Correction audit](results/accuracy-review/capsid_vertex_audit.jsonl)
- [Naming equivalents](naming_equivalences.tsv)
- [External-panel selection](external_panel/selection.json)
- [Checksums](results/accuracy-review/SHA256SUMS)
