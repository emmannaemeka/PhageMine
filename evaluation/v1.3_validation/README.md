# Fresh annotation validation

[Fresh seven-genome workflow](https://github.com/emmannaemeka/PhageMine/actions/runs/37163696543)
completed all seven jobs successfully. PhageMine, Pharokka 1.10.1 and Prokka
1.14.6 each generated new annotations from the same nucleotide-only inputs.
PhageMine and Pharokka shared PHROGs v4 from the Pharokka 1.8.0 distribution.
PhageMine used PHROGs only; Prokka used its packaged viral databases. RNA
features were disabled for this CDS comparison. Commands and database/package
versions are retained in the workflow artifacts.

The panel contains 462,919 nucleotides and 700 simple reference CDS. Eleven
compound reference CDS are excluded. It is a familiar diagnostic panel, not an
independent held-out test. Product names in the reference files never enter the
FASTA annotation inputs.

## Findings and fixes

Fresh evidence exposed partial-profile naming and correlated-profile voting
failures. Two short T4 proteins were assigned the name of a whole enzyme from
matches covering only 7.1/8.7% of its profile. Three partial tail-fiber profiles
outvoted a stronger full-length connector match. Fusion rules 1.10 withhold the
former names and select the supported connector description.

The code audit also identified gene/EC qualifiers borrowed from rejected or
unselected product hypotheses. Qualifiers now require strong eligible evidence
supporting the selected whole-protein product. Missing annotations and
unknown-function family names cannot establish a named function. Synthetic
regression tests exercise qualifier transfer; this PHROGs-only benchmark does
not independently validate Swiss-Prot gene-symbol transfer.

Current rules were applied to the newly generated raw evidence using
`phagemine reclassify`. This is explicitly a decision audit, not another database
search. All seven corrected exports passed artifact validation. Raw evidence
hashes and gene coordinates remained unchanged. **41 product assignments
changed**; these are review candidates, not 41 independently proven corrections.

| Decision set | Predicted CDS | Exact simple reference models | Named assertions at named reference loci | Literal normalized product agreements | Unresolved named differences |
| --- | ---: | ---: | ---: | ---: | ---: |
| PhageMine fresh, rules 1.9 | 804 | 588 | 346 | 279 | 67 |
| Pharokka fresh | 804 | 588 | 341 | 303 | 38 |
| Prokka fresh | 678 | 616 | 268 | 35 | 233 |
| PhageMine fresh evidence, rules 1.10 | 804 | 588 | 335 | 276 | 59 |

Counts use 700 simple reference models, including 465 named reference loci.
Some named reference loci lack an exact predicted model. Literal name agreement
is not functional accuracy: synonyms and specificity differences require review.
Functional precision/recall and superiority claims remain withheld. Search
policies and Prokka's evidence databases differ, so this is not a controlled
comparison of adjudication alone.

The partial-profile rule increases abstention at exact named reference models
from 44 to 55. It may withhold valid names from genuine protein fragments.
That tradeoff requires independent review, rather than optimizing name counts
against these seven familiar genomes.

A capsid vertex vs major-capsid distinction remains unresolved: PhageMine and
Pharokka both inherit a coarse PHROG family assignment at the T4 locus. The
raw evidence for this case is included. Copying a comparator's label would not
resolve it. More discriminating curated evidence and an independent held-out
panel remain necessary before claiming better overall annotation accuracy.
Phold was not run. No wet-laboratory validation is required for this software
audit or claimed by these results.

## Reviewable data

The [paired accuracy review](ACCURACY_REVIEW.md) separates differences between
tools from differences against the reference and records the review protocol.

- [Manifest and artifact hashes](results/fresh-seven/manifest.json)
- [Complete summary](results/fresh-seven/summary.tsv)
- [41 changed product assignments](results/fresh-seven/decision_changes.tsv)
- [Raw evidence for changed assignments and the capsid-vertex case](results/fresh-seven/review_evidence.json)
- [Fresh predictions, corrected predictions and reference tables](results/fresh-seven/inputs/)

The initial serial run's completed T4 results are separately documented in
`results/fresh-t4-diagnostic.json`; its remaining batch was canceled when the
panel moved to parallel jobs. They are not substituted for the completed
seven-genome panel.

## Reproduction

Use `.github/workflows/fresh-benchmark.yml` to generate new annotations. To
audit current decision rules, extract the seven workflow artifacts and run:

```bash
python evaluation/v1.3_validation/scripts/evaluate_fresh_results.py \
  --results /path/to/NC_000866.4 /path/to/NC_000929.1 \
  /path/to/NC_001416.1 /path/to/NC_001422.1 /path/to/NC_001604.1 \
  /path/to/NC_002371.2 /path/to/NC_005859.1 \
  --output fresh-decision-audit
```

The audit rejects mixed databases/source commits, archived-run declarations,
changed scientific artifacts, duplicate genomes and changed raw evidence or
coordinates. GitHub omitted hidden Matplotlib font caches from the first
artifacts; those disposable rendering-cache omissions are explicitly recorded.
Scientific artifacts must pass their recorded checksums.

`scripts/run_archived.py` remains a historical scoring diagnostic. Its results
are not fresh annotation measurements and are not used for the table above.
