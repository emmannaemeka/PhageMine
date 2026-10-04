# Annotation naming policy

Gene symbols (for example `polA`), locus identifiers (`PM_000001`), and protein
products (`DNA polymerase`) are different fields. PhageMine does not infer a
gene symbol from a product name or copy the name of a neighboring gene.

The naming audit inspected Prokka's `cleanup_product` and qualifier transfer,
and Pharokka 1.10.1's PHROG annotation-table merge:

- [Prokka naming implementation](https://github.com/tseemann/prokka/blob/master/bin/prokka)
- [Pharokka 1.10.1 annotation implementation](https://github.com/gbouras13/pharokka/blob/v1.10.1/src/pharokka/post_processing.py)

The shared principles are evidence-backed product transfer, explicit separation
of identifiers and functions, and hypothetical products when no function is
established. Pharokka's internal `gene` column also identifies query proteins;
it is not automatically a biological gene symbol. Neither tool's output is
treated as ground truth merely because the tools disagree.

Fusion rules 1.11 implement these naming safeguards:

1. Missing descriptions (`NA`, `N/A`, `nan`, `-`) and unknown-function families
   (`DUF`, `UPF`) cannot establish a named function. Raw descriptions remain in
   evidence; conservation may still be reported.
2. Gene and EC qualifiers require a strong or experimental accepted record
   with the same normalized function as the selected product. Domain-only
   records, rejected products, conflicts with no selected product, and
   member-specific identifiers from a broader family rewrite cannot transfer
   their qualifiers.
3. Conflicting explicit identifiers are withheld. No new symbol is invented.
4. Informative descriptions containing `conserved` are preserved. Uncertainty
   words such as `putative` are preserved. Prokka's broad cleanup rules are not
   copied wholesale because they can discard useful phage-specific detail.
5. A PHROGs/VOGDB match covering less than half of the reference profile cannot
   supply a whole-protein product name. Its raw alignment and conservation
   evidence remain available. This is a conservative rule, not an empirically
   calibrated probability, and may withhold names from genuine fragments.
6. Multiple profiles from one database count as one source for consensus.
   Repeating a profile cannot outvote a stronger distinct function. Score
   margins compare different product hypotheses, not duplicate records of the
   winning hypothesis. Different sources are not assumed statistically independent.
7. A reviewed whole-protein name keeps its specificity and uncertainty when
   a generic domain rule is also present. Domain rules cannot choose between
   contradictory reviewed whole-protein records; those cases remain unresolved.

Use the most informative role supported by the selected evidence. Broader family
labels cannot establish a specific family member, and unknown proteins must not
be named by copying a neighbouring gene or a comparator's assertion. Explicit
curation preserves the original automatic result and records its evidence.

The reviewed T4 sequence matching UniProt P19896 supports **capsid vertex
protein**, rather than major capsid protein. This user-approved correction is
retained in an audited curation file and reported separately from automatic
benchmark performance. It is not a locus-specific rule in the annotation engine.

Fresh T4 evidence exposed two specific failures of the preceding rules:
61/67-aa proteins were assigned a whole enzyme name from matches covering only
7.1/8.7% of its profile; and three partial tail-fiber profiles outvoted a stronger
full-length connector profile. Rules 1.10 withhold the former names and select
the supported connector description. This diagnoses decision-rule failures;
it does not establish general annotation accuracy.

## Fresh comparison

The `Fresh annotation comparison` GitHub Actions workflow runs current
PhageMine, Pharokka 1.10.1 and Prokka 1.14.6 on nucleotide-only FASTA inputs
derived from the seven reference genomes. All predictions are generated anew;
archived predictions are not annotation inputs.

PhageMine and Pharokka share PHROGs v4 distributed with Pharokka database 1.8.0.
PhageMine uses both PHROGs backends; Pharokka uses its standard sequence/profile
searches. Prokka uses its packaged viral databases. Search policies differ and
are retained in command logs and adapter manifests. This tests a PHROGs-only
PhageMine configuration, not the additional Swiss-Prot, Pfam and VOGDB evidence.

Artifacts include input/output checksums, package versions, raw runs, exact-locus
product comparisons, and `gene_name_review.tsv`. Gene-symbol differences include
missing symbols and locus-specific aliases; they are review candidates, not
counts of proven incorrect names. Compound reference CDS are excluded from the
contiguous-coordinate analysis. Failed searches stop the benchmark rather than
being scored as missing functions.

The seven familiar genomes are a diagnostic panel, not a held-out validation
set. Fresh predictions do not make reference annotations independent or prove
accuracy superiority. Unresolved naming differences need evidence-backed review.

To audit updated decision rules without repeating identical searches, extract
the per-genome fresh artifacts and run:

```bash
python evaluation/v1.3_validation/scripts/evaluate_fresh_results.py \
  --results /path/to/extracted/genome1 /path/to/extracted/genome2 \
  --output fresh-decision-audit
```

This verifies completed fresh-run status, input/artifact checksums, matching
database hashes and source commits. It regenerates PhageMine products using
current rules, validates exports, preserves raw evidence and coordinates, and
compares the original fresh decisions, corrected decisions, Pharokka and Prokka.
`decision_changes.tsv` records changes for review. Corrected decisions are
explicitly labeled as reclassification of fresh evidence, not another search run.
