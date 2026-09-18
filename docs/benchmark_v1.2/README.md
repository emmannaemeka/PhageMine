# Public benchmark package

This compact package contains the frozen seven-phage benchmark summaries and
publication figures. It deliberately excludes raw tool outputs, annotation
logs, and third-party databases. The complete local evidence remains preserved
outside the public source package. The small subset of frozen reference and
annotation files needed for functional validation is now archived in
`evaluation/v1.2_functional/`.

Tables use one-based inclusive coordinates and explicit strict/relaxed
structural definitions. Reference provenance records the accession-specific
GenBank inputs and checksums. Figures are derived from the TSV tables.

Verify the package with `shasum -a 256 -c docs/benchmark_v1.2/SHA256SUMS` from the repository root.

## Functional-yield correction (18 September 2026)

The previous per-genome table incorrectly reported zero named PhageMine
products for all seven genomes. The historical structural generator read
product labels from `genes.gff3`, which did not contain PhageMine's functional
products, instead of joining `annotation.tsv`. Those zeros were invalid; the
validated aggregate of **439** was correct and is unchanged.

Authentic frozen annotation tables were recovered from the original corrected
benchmark workspace. PhageMine counts are T4 **144**, Lambda **56**, T7 **50**,
T5 **95**, PhiX174 **8**, P22 **48**, Mu **38** (sum **439**). Pharokka and Prokka
per-genome counts were independently verified against their frozen GFFs and
are unchanged. The per-genome table now records source paths and SHA256 hashes.
The compact source files and recovery provenance are archived in
`evaluation/v1.2_functional/`; no annotation tool was rerun.

The functional-yield plot is regenerated in PNG, SVG and PDF directly from
`tables/functional_yield_totals.tsv`, with assertions for 439, 398 and 328.
Structural results and structural figures are unchanged. The figure script
regenerates functional yield only; it no longer emits an unrelated structural
plot using a different aggregate averaging rule.

Functional yield does not establish accuracy. See the new functional benchmark
for explicit reference concordance, unresolved manual review and limitations.
