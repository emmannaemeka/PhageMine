# Frozen v1.2 functional reference-concordance benchmark

Status: **automated analysis complete; biological adjudication incomplete**.
No tool was rerun. Missing judgments are `UNRESOLVED_REVIEW_REQUIRED`, and
metrics that depend on them are `NA` (JSON `null`), never invented zeros.

## Reproduction

From the repository root, use Python 3.10+, Biopython 1.85 and Matplotlib 3.x:

```sh
python evaluation/v1.2_functional/scripts/recover_yield.py
python evaluation/v1.2_functional/scripts/run_benchmark.py
python evaluation/v1.2_functional/scripts/make_figures.py
python docs/benchmark_v1.2/scripts/make_figures.py
python evaluation/v1.2_functional/scripts/checksums.py
shasum -a 256 -c evaluation/v1.2_functional/SHA256SUMS
shasum -a 256 -c docs/benchmark_v1.2/SHA256SUMS
```

The scorer imports `src/phagemine/functional_benchmark.py` directly. It does
not use the older `benchmark.py` heuristic product-comparison API, which
contains broader terminology substitutions and treats unresolved differences
as disagreements. That API is not suitable as biological ground truth here.
No network, databases or annotation executables are needed for this analysis.
Use `--output /path/to/new/results` to verify regeneration without modifying
this snapshot. Input SHA256 assertions and yield assertions fail on drift.

## Evidence and the zero-count correction

The original local workspace was
`/Users/emmanuelnnadi/Documents/PhageMine_v1.2_development/evaluation/v1.2_tool_comparison/`.
The source inventory and individual hashes are in `source_provenance.tsv`.
The archived structural generator and its method provenance are retained in
`provenance/` as historical evidence, not as scripts to execute.

* References: original `inputs/<accession>.gb`, copied byte-for-byte. All seven
  checksums equal the published reference provenance. Their original retrieval
  date was not recorded; the archive date is not presented as retrieval date.
* PhageMine: `correction_validation_20260916/outputs/phagemine/<accession>/`
  containing `annotation.tsv`, `genes.gff3`, `functional_classification.tsv`
  and `run_manifest.json`.
* Pharokka: `raw_outputs/pharokka/<accession>/<accession>.gff`.
* Prokka, used only to verify Task A:
  `correction_validation_20260916/outputs/prokka/<accession>/<accession>.gff`.

The previous structural script read PhageMine product names from the GFF
instead of the already-loaded annotation table. The GFF has coordinates but
no functional products. Correctly joining the annotation table recovers
144/56/50/95/8/48/38 named calls, totaling 439; Pharokka 398 and Prokka 328 are
also reproduced. Each annotation join asserts ID and coordinates. The
original lexical yield definition is preserved separately from the stricter
reference-informativeness rule. No aggregate count is expanded into synthetic
individual records. Source tool evidence is retained for review, but agreement
among annotation databases is not assumed to be independent validation.

## Reference truth and matching

All 711 reference CDS features are retained, with full GenBank location strings
and qualifiers in `functional_truth_by_locus.tsv`. Informative reference labels
exclude hypothetical, uncharacterized/uncharacterised, unknown and unassigned
labels, including decorated versions; generic “protein” and “phage protein”
labels are also excluded. Gene symbols alone are not silently expanded into
functions. This is reference concordance, not experimental validation.

Primary scoring uses accession version, one-based inclusive start/end and
strand. Both tools have 804 identical coordinate tuples. There are 588 exact
matches, including 390 evaluable named-reference loci and 198 non-evaluable
reference loci. Metrics are conditional on the matched evaluable set; they
are not genome-wide annotation recall. Unmatched CDSs are reported separately.

Secondary scoring locks exact matches first and applies the published
same-strand overlap thresholds: >=50% of reference span and >=20% of prediction
span. Only unique one-to-one residual links are accepted. Ambiguous or
many-to-one overlaps are recorded in `matching_audit.tsv`, not forced into
functional matches. Thus this functional subset does not recreate the
structural relaxed TP count. It has 652 matches, 435 evaluable.

Eleven reference CDSs have compound locations (including circular-origin
joins). Their full locations are retained, but they are excluded from functional
matching: the historical structural parser's bounding intervals can span
unrelated CDSs. This exclusion leaves the 588 primary exact matches unchanged.
It limits the secondary analysis. No structural table or score is revised.

## Classification and blinded review

Automatic decisions are restricted to normalized full-product equality,
explicit full-term synonym pairs in `synonym_rules.tsv`, abstention, and
non-evaluable reference. Normalization changes case/punctuation/spacing only;
it does not delete putative/probable qualifiers or merge subunits. The sole
initial synonym rule is the protocol's major-head/major-capsid equivalence.
No fuzzy similarity or LLM judgment is used as truth.

All other semantic cases await review. The eight final categories are
`EXACT_PRODUCT_AGREEMENT`, `EQUIVALENT_FUNCTION`, `COMPATIBLE_BUT_BROADER`,
`COMPATIBLE_BUT_MORE_SPECIFIC`, `ABSTENTION`, `UNSUPPORTED_SPECIFICITY`,
`GENUINE_FUNCTIONAL_DISAGREEMENT`, `NON_EVALUABLE_REFERENCE`.
`UNRESOLVED_REVIEW_REQUIRED` is a workflow state, not a ninth biological outcome.

Give reviewers **only** `manual_adjudication_blinded.tsv`, not the key, the
prediction archives, or the unblinded results. A/B identity is independently
shuffled per locus/matching mode with a fixed seed. The mapping is separately
stored in `manual_adjudication_key.tsv` for reproducibility. Public availability
of that key means blinding requires procedural separation; it is not secrecy.
The 621 review rows include exact and relaxed strata separately. Duplicate
loci across strata should receive consistent biological judgments.

Copy the review TSV, fill only unresolved sides, use the exact category names,
and provide a rationale in `reviewer_notes`. More-specific assignments require
an `evidence:` citation documenting independent support; the software checks
that a citation is supplied, not its biological validity. A qualified reviewer
must assess its independence. Broader calls must preserve the reference
function; unsupported specificity is distinguished from an actual biological
conflict. Reviewers must not assume an extra named call is correct.

```sh
python evaluation/v1.2_functional/scripts/run_benchmark.py \
  --adjudications /path/to/completed_review.tsv --output /path/to/reviewed_results
```

Review ingestion validates IDs, immutable product/coordinate fields, categories
and justifications. It cannot overwrite deterministic decisions. No review is
silently applied; the submitted file checksum is recorded. Keep the completed
review file with any subsequently published adjudicated run.

## Metrics and statistics

Denominators are matched evaluable named-reference loci (`N`) and their named
assertions (`A`). Non-evaluable references are never counted as correct.

* Coverage = A/N; abstention = abstentions/N.
* Strict precision = (exact + equivalent)/A.
* Expanded precision = (exact + equivalent + broader + justified more-specific)/A.
* Functional recall = (exact + equivalent)/N; F1 is their harmonic mean.
* Unsupported specificity and disagreement rates use A as denominator.

`functional_summary.tsv` contains category counts, numerator/denominator,
proportions, percentages and Wilson 95% intervals where applicable. With
unresolved named assertions, correctness/error metrics and their confidence
intervals are withheld. Identification bounds show the range attainable if
none/all unresolved judgments enter that metric; these are **not confidence
intervals**. Zeros in unreviewed semantic-category counts mean no judgments
entered, not demonstrated absence of those errors. Undefined denominators
produce NA. F1 confidence intervals are not computed.

The paired exact McNemar calculation is implemented for adjudicated correct /
not-correct outcomes on the common evaluable loci, but no p-value is reported
while any outcome is unresolved. Even after review, loci within a phage are
not independent: locus-binomial intervals and the unclustered paired test are
descriptive, and cannot support population-wide claims. Seven genomes do not
provide a strong basis for generalization or a forced significance claim.

## Extra assignments and limits

`phagemine_only_functional_assignments.tsv` contains **75**, and the reverse
file **34**, observed exclusive named assignments. Their difference is 41;
neither exclusive set has size 41. Among PhageMine's 75: 4 exact product
agreements, 32 pending semantic review, 21 non-evaluable references and 18
without an exact reference match. The reference-supported proportion among
36 evaluable exact loci and its CI await review. Its current identification
range is 4/36 to 36/36; it is not a point estimate. All-75 validation is limited
by the other 39 loci, even if review of the 32 is completed.

These classic phages may be represented directly or indirectly in PHROGs,
VOGDB and Swiss-Prot. Agreement measures reference-concordant annotation, not
independent performance on novel phages. No universal superiority claim is
permitted. A second stage should use temporal or sequence-held-out genomes,
with database snapshot dates and sequence exclusion verified. A local temporal
panel exists at `evaluation/v1.2_finalization/temporal_primary_comparison.tsv`
in the original development workspace; its accession list is archived as
`provenance/temporal_candidates.tsv`. These are candidates, not independently
validated held-out genomes; they are excluded from every result here.
