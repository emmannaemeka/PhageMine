# Full-system naming comparison with Prokka

PhageMine now recovers **20 of the 25 shared genes** where Prokka supplied a
name and the previous full-system adjudication withheld one. All 20 recovered
names are consistent with their exact reviewed sequence records. Five remain
withheld. Gene calling, all 978 gene coordinates and their protein sequences
are unchanged. All ten regenerated outputs pass integrity checks.

This comparison uses PhageMine's final products from PHROGs, Swiss-Prot, Pfam
and VOGDB evidence, with PMFDB/INPHARED comparative analysis completed. It uses
729 exactly shared loci and the earlier Prokka 1.14.6 outputs from identical
nucleotide inputs. Database versions, parameters and workloads differ; this
is a diagnostic comparison, not a controlled accuracy or speed trial.

## What changed

Rules 1.14 prefer a complete, ungapped 100%-identity reviewed phage alignment
over differently named, less-identical homologues. The complete query,
reference and alignment lengths must agree. This alignment-based priority is
not itself an independently verified sequence match or experimental function.
Distinct informative names from equally complete full-identity reviewed
records still require review. Less-identical alternatives remain visible, and
uncertainty words are retained. Perfect partial matches receive no such priority.

This resolves unnecessary abstentions for portal, terminase, capsid decoration,
internal virion and other proteins. It also distinguishes the reviewed T7 major
capsid sequence from a related minor-capsid isoform. Unknown-function records
cannot veto an informative candidate merely because they are reviewed.

Numeric gene wrappers, prophage-origin text attached to an identifier, and
names consisting only of a database-entry symbol cannot establish function.
Informative enzyme/structural roles remain eligible. Domain-only evidence stays
separate from whole-protein naming. There are no genome/locus-specific engine
lookup rules or comparator-label transfers.

The final full-system output retains **capsid vertex protein** automatically
for the T4 sequence matching reviewed P19896.

## Evidence review

The reviewed cohort is the union of 273 initial naming differences and 253
current differences: **274 loci**, including cases resolved by the new rules.
This captures every current naming difference at shared coordinates. Unshared
gene models are reported separately and are not scored as unnamed proteins.

The frozen review contains 151 unique versioned reviewed records. Protein
sequence and record checksums bind every supported decision to its exact
record. Relevant record fields are frozen without rewriting; full source-record
checksums are retained in the projection manifest. There are **154 exact reviewed sequence cases** and **120 without an
exact reviewed match**. High identity is not silently treated as an exact
sequence match. Ten additional candidate records were fetched from UniProt
for this review; candidate accessions alone did not establish sequence identity.

| Review outcome | PhageMine | Prokka |
| --- | ---: | ---: |
| Named role consistent with reviewed record | 111 | 45 |
| Compatible broad role | 2 | 6 |
| Specificity not established | 4 | 1 |
| Assertion not established by record | 32 | 2 |
| Function/name conflicts with record | 0 | 2 |
| No functional role supplied | 5 | 98 |
| No exact reviewed sequence match | 120 | 120 |

These counts are **not precision estimates or evidence of overall superiority**.
The sample selects differences, including many PhageMine-named/Prokka-unnamed
cases. Direct agreement with an informative recommended/alternative record name
supports nomenclature; other decisions explicitly cite FUNCTION comments,
record identifiers, isoforms or keywords. The reviewer is Codex, not an
independent domain expert, and the original papers were not independently
reviewed. Annotation and review sources overlap. Reviewed records can contain
inferred functions. Missing support does not prove an assertion false.

Concrete findings:

- T4 P04527 is the small clamp-loader subunit, gene 62. Prokka calls the exact
  sequence a gp44 subunit; PhageMine keeps the reviewed small-subunit name.
- T7 P19726 is the displayed major capsid isoform. The record distinguishes it
  from the external minor isoform and their 90/10 ratio (PubMed:20962334).
  PhageMine selects the major role; Prokka assigns the minor role.
- P22 P26746 is both a head-to-tail adapter and a peptidoglycan hydrolase.
  Both tools' descriptions are supported by the record; differing wording is
  not an error (PubMed:16970964,14763988).
- Putative/probable and member-number differences were reviewed case by case.
  No broad synonym rule was added to inflate the benchmark.

The five remaining Prokka-named/PhageMine-unnamed cases include partial enzyme
matches and two sequences whose exact reviewed entries are uncharacterized.
Those specific functional claims remain unresolved rather than being copied
from Prokka. All 32 PhageMine assertions not established by their record and
four specificity cases also remain review candidates. In particular, absence
of a FUNCTION comment must not automatically veto other accepted evidence.

## Reproduce the review

Extract the full-system artifact from
[run 37186254566](https://github.com/emmannaemeka/PhageMine/actions/runs/37186254566).
The searches were generated fresh; the naming-rule audit reuses those validated
raw alignments without changing genes or sequences.

```bash
python evaluation/v1.3_validation/scripts/review_full_system.py \
  --source full-results --output reviewed-full-results

python evaluation/v1.3_validation/scripts/complete_accuracy_review.py \
  --paired evaluation/v1.3_validation/results/prokka-review/selected_cases.tsv \
  --runs-root reviewed-full-results/runs \
  --snapshot evaluation/v1.3_validation/results/prokka-review/reviewed_records.json \
  --decisions evaluation/v1.3_validation/results/prokka-review/decisions.tsv \
  --right-tool Prokka --output reviewed-prokka-cases
```

The second command rejects changed names, sequences, records, missing decisions
and duplicate loci. It requires the rule version/source checksums recorded with
this review; future changed outputs need new decisions.

- [All 20 recovered names and their records](results/prokka-review/recovered_names.tsv)
- [Comparison provenance](results/prokka-review/comparison_manifest.json)
- [Every selected case and decision](results/prokka-review/reviewed_cases.tsv)
- [Review interpretation](results/prokka-review/review_summary.json)
- [Initial paired output](results/prokka-review/paired_initial.tsv)
- [Current paired output](results/prokka-review/paired_current.tsv)
- [Record projection and source hashes](results/prokka-review/record_projection_manifest.json)
- [Versioned reviewed records](results/prokka-review/reviewed_records.json)
- [Review checksums](results/prokka-review/SHA256SUMS)
- [Full-system outputs and reference counts](BENCHMARK_REVIEW.md)
