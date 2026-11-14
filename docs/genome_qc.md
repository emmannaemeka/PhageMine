# Genome QC and structural curation

PhageMine implements a computational checklist aligned with the 14 questions in
Turner et al., *Phage Annotation Guide* (2021),
https://pmc.ncbi.nlm.nih.gov/articles/PMC8785237/.
A PASS is scoped to the described computational check. The report never asserts
physical genome completeness, therapeutic suitability or INSDC acceptance.

## Pipeline integration

PHANOTATE uses an explicit 90-nt minimum. Pyrodigal providers use 90-nt normal
and edge minima. After same-stop Prodigal-GV selection, `start_refinement.py`
enumerates table-11 starts up to the preceding in-frame stop, scores SD-like
motifs, spacing, codon preference and optional full-length DIAMOND N-terminal
support. It preserves gene IDs and records every adjustment in
`gene_calls/start_refinement/changes.json`. All subsequent functional evidence
is searched against the refined proteins. No reference coordinates enter this
runtime path. Existing explicitly partial genes are not forced to become complete.

Rare CTG/ATT/ATC/ATA starts are permitted and flagged for review. SD scoring is a
heuristic, not evidence that all phages require an SD motif. An extension without
an SD motif can be supported by compact spacing and is labelled separately for
review. Candidate scores and unchanged alternatives remain available in JSON.

The short-ORF scanner evaluates 75–300-nt candidates on both strands; lengths
include the stop codon (75 nt encodes 24 amino acids plus the stop). Rescue requires
initiation or homology support **and** an unannotated gap. At most four bases may
overlap either neighboring CDS boundary; nested and antisense overlapping genes
are never automatically rescued. Supported rare starts are eligible. Candidates
rejected for overlap remain visible for manual curation. This deliberately cannot
recover every genuine overlapping phage gene.

The default heuristic was evaluated on PZ583544.1. Improvement on this benchmark
is not an independent estimate of accuracy on other phages.

## Reports and commands

Normal export writes `phagemine_qc_report.json` and `.md` in the run directory
and submission package. The report reads final GenBank features and checks
ambiguous bases, homopolymers (>=10 bases), exact terminal-repeat candidates
(20–10,000 nt), internal copies of detected terminal motifs, CDS codons/frame,
RNA/CDS overlaps, naming and evidence diversity. Repeat detection is an exact
end-repeat screen, not an exhaustive internal repeat or assembly-completeness
analysis. Ordered adjacent homology pieces in different frames are flagged only
when subject alignment coordinates are present; domain names alone are inadequate.

```sh
phagemine genome-qc RESULTS_DIR
phagemine genome-qc RESULTS_DIR --safety-db amr=/path/amr.dmnd \
  --safety-db toxin_virulence=/path/virulence.dmnd \
  --safety-db lysogeny=/path/lysogeny.dmnd --threads 8
```

Dedicated screens record database/query SHA-256, commands and thresholds. Default
review thresholds are 40% identity, 80% coverage of both proteins, E <= 1e-5.
They do not substitute for AMRFinderPlus/CARD-specific curated models or laboratory
safety testing. A screen is invalidated if its protein input changes. Missing
categories, failed searches and absent physical/read evidence remain explicit.

Existing tRNAscan-SE annotation is retained. ARAGORN, if installed alongside it,
adds a linear, table-11 tmRNA scan. Unrecognized/permuted RNA locations are not
silently emitted as simple intervals. RNA overlaps are review findings, not an
automatic deletion rule. Computational regulatory-feature inputs can be checked
with `curate_regulatory`; CDS-overlapping predictions are rejected. PhageMine does
not invent promoters or terminators when a predictor is not configured.

```sh
phagemine assemble --reads1 R1.fastq.gz --reads2 R2.fastq.gz \
  --assembler shovill --output new_assembly --threads 8
phagemine run new_assembly/contigs.fasta --output annotation
```

The optional assembly command also accepts `--assembler spades`. Shovill uses
`--depth 100`; that is a subsampling cap, not a verified 25–100x coverage result.
Assemblies never overwrite an existing output directory. Raw-read hashes,
commands and logs are retained. No raw reads were supplied for the PZ benchmark.

## Benchmark

```sh
python -m phagemine.cds_benchmark --prediction RESULTS/genbank_submission/genome.gbk \
  --reference PZ583544.1.gb --output benchmark.json
# GFF3 also accepts --fasta sequence.fasta
```

Matching is one-to-one by contig, strand and stop, prioritizing exact boundaries.
The report lists exact matches, stop-only matches with codons, missed reference
features and unmatched predictions (benchmark false positives). Reference data
belong only in evaluation inputs and tests. Existing IDs are never renumbered to
improve benchmark metrics. Local `PM_000001` syntax does not establish that its
prefix is registered for an INSDC submission.
