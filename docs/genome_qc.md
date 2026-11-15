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

## Spanin subtype normalization

After functional evidence fusion, `spanins.py` normalizes `rz spanin` and
`rz-like spanin` to `Rz-like spanin`. A same-contig, same-strand two-component
cassette can promote the downstream spanin to `Rz1-like spanin` when it overlaps
in the +1 relative reading frame, or is separated by at most 15 nt. It must also
have the requested `[LVI][ASTVI][GAS]C` motif entirely within its first 30 residues,
or an accepted Rz1-family match. The initial explicit family policy includes
PHROGs `phrog_5672`; explicit Rz1/o-spanin/outer-membrane spanin descriptions from
other evidence sources are also recognized. Family support requires REAL,
supporting evidence, E <= 1e-5 and >=50% query coverage. These are subtype-evidence
requirements, not new gene-calling thresholds. Competing cassette partners are
left unresolved. Unimolecular spanins and wrong-frame overlaps are not forced
into the two-component model.

The decision and evidence are retained in the classification record and protein
provenance. Coordinates, sequences and caller thresholds do not change. This is
a computational product-name assignment, not proof of membrane localization.
The same rule works in transcription direction on the negative strand. Some
biological spanin pairs have other frame arrangements; the narrow requested
rule deliberately leaves those unresolved.

The six-genome evaluation driver is `evaluation/multi_phage/benchmark.py`.
Fixtures are external NCBI GenBank/FASTA files, and runtime input is FASTA only.
The benchmark retains ordered joined intervals and does not credit a simple CDS
spanning an intron or frameshift as an exact structural match. It reserves exact
matches before assigning other predictions sharing a stop. All missed features,
extra predictions, tRNA matches and strict codon failures remain in its JSON.
