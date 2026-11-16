# Benchmark-driven structural QC repairs

These rules use sequence, caller agreement and protein evidence. Reference coordinates are used only by evaluation scripts.

* Partial boundaries within two bases of a contig edge extend to the physical edge in biological strand orientation. Added 5′ bases are skipped using `codon_start`; incomplete 3′ codons are excluded from translation. The same coordinates and frame offset are emitted in GenBank, feature tables and GFF3. This does not convert internal truncations into contig-edge features.
* Translation uses Biopython table 11, including IUPAC resolution (for example GTK → V). Reverse complements also preserve full IUPAC semantics. Only one terminal stop is removed, so internal stops remain detectable.
* tRNAscan-SE and optional ARAGORN execute before short-ORF rescue. Their validated RNA intervals exclude candidates on either strand. RNA records are attached after the output files exist, without running the searches twice.
* Ordinary rescue remains restricted to intergenic gaps. A separate post-evidence pass examines same-strand downstream +1-frame overlap candidates of 75–300 nt. RNA overlap always blocks promotion. Rz1 recovery requires an annotated inner-spanin context, a strong RBS (score ≥8), and a lipobox or Rz1-family support. Other overlap promotions require real, accepted PHROGs/VOGDB evidence with E ≤1e-5 and query coverage ≥70%, from a family distinct from all overlapping proteins. Antisense and other-frame candidates remain in review. Candidate evidence and decisions are recorded in `overlap_rescue.json`.
* Above 58% GC, a smoothstep adjustment over 58–70% GC reduces the ATG-to-GTG/TTG prior gap by up to 80%. Moderate-GC priors are unchanged. SD prevalence is measured only at complete starts where the two callers originally agree; at least ten anchors are required. In high-GC genomes with <35% strong-SD anchors, SD scoring is disabled, avoiding selection against otherwise supported SD-free starts. An absent SD motif does not establish leaderless transcription.

## Evaluation

Run `scripts/verify_codons.py genome.gbk --tbl features.tbl --extended-starts --require-table 11` for independent exported-feature verification. Use `evaluation/multi_phage/benchmark.py` to compare against externally supplied GenBank references. Inspect exact, shared-stop, missed and extra predictions separately; a legal CDS does not prove a biologically correct gene or start.

The baseline Bxb1 audit contains nine ATG, four GTG and two TTG reference starts among 15 shared-stop mismatches. Thirteen have an SD score of zero, one scores 4.8 and one scores 8. However, 36 of 61 caller-agreement anchors have strong SD motifs, so Bxb1 does not activate the <35% rule. The data do not justify relaxing that threshold to force benchmark matches.

Codon verification and local unit tests do not establish NCBI submission acceptance. Required source/submitter metadata and official validation must be assessed separately.
