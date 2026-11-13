# Annotation reliability review

Fresh single-record DNA annotation runs now write `annotation_reliability/index.html`,
linked from the main report, plus machine-readable JSON and TSV files.
The checks are observational: no final CDS coordinates, products, translations,
or whole-genome similarity methods are replaced.

## Automatic checks

- Coding capacity is the union of valid, inclusive CDS intervals divided by genome
  length. Overlaps are counted once. Values below 85% produce a review warning in
  `quality_control.json`. This pragmatic threshold is not a biological cutoff:
  low coverage may reflect assembly errors, incomplete input, missed genes,
  noncoding sequence, or alternative coding. It does not establish recoding.
- Successful native Pyrodigal-gv observations expose translation tables and
  secondary coding capacity. Non-table-11 predictions are flagged for review.
  Missing or failed secondary observations are `NOT_ASSESSED`; stale output from
  a failed secondary run is not used. This is native-model inference, not a new
  exhaustive test of genetic codes. Final translations are unchanged.
- CDSs shorter than **100 nucleotides**, not amino acids, receive a review table
  with caller agreement and available functional evidence sources. Secondary-only
  short calls are flagged as possible omissions or disputed boundaries. These are
  candidates, not automatically accepted genes. Genes missed by both callers
  cannot be detected by this comparison; native minimum-length limits still apply.

## Annotated reference checks

Supply annotated reference genomes obtained from RefSeq or an appropriate,
independently verified ICTV reference set:

```bash
phagemine run genome.fasta --output results \
  --reference-genbank references.gb \
  --reference-min-coverage 70 --reference-min-identity 80
```

Repeat `--reference-genbank` for more files. BLASTN is required for this optional
step; `--reference-blastn` selects an executable. The existing INPHARED workflow
continues independently. This feature does not download GenBank annotations or
verify ICTV exemplar status automatically. RefSeq accession prefixes are labelled
as provenance clues; all ICTV statuses remain `NOT_VERIFIED`.

The screening identity is explicitly defined as identical A/C/G/T positions over
aligned query bases, selecting the highest-scoring HSP for overlapping query
positions. Query bases opposite subject gaps count in the denominator; subject
insertions do not. This is not BLAST's per-HSP `pident` or whole-genome similarity.
All references, including those failing the adjustable 70%/80% screen, appear in
`reference_candidates.tsv` and a coverage-versus-identity SVG. Thresholds are
screening heuristics; the plot does not claim statistical outlier significance.

For passing references, simple CDS features with exact boundaries are projected
through a single alignment. Exact boundary agreement, disagreement, possible
missing CDSs, ambiguous mappings and indels are reported separately. Compound,
partial, split-HSP and unaligned features are reported as unmappable. Reference
support does not prove gene function and no annotations are transferred.

The reference report records input hashes and commands. Duplicate reference IDs,
invalid input, BLAST failures and timeouts are visible as a failed optional review;
annotation remains usable. Missing BLAST is reported as unavailable.

This implementation adds no VIRIDIC dependency and changes no taxonomy thresholds.
Review remains restricted to single-record DNA genomes. Reference arguments on
RNA or segmented runs are rejected explicitly.
