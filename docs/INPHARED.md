# INPHARED comparative genomic context

PhageMine can add nearest-reference phage context after annotation. The
workflow is:

```text
input genome → PHANOTATE CDSs → PHROGs/VOGDB/Pfam/Swiss-Prot evidence
→ annotated genome → INPHARED comparative genomic context
→ nearest/reference relationships and supported metadata
```

## Install and run

```bash
phagemine databases install inphared
phagemine doctor
phagemine run genome.fasta --output results/genome
```

For a cohort-level comparison:

```bash
phagemine batch genomes/ --output results/cohort --mode discover --evidence core
```

The pinned resource is INPHARED release **2026-04-07**. Installation downloads
the genome FASTA and metadata table, verifies provider MD5 checksums, creates a
per-genome Mash sketch with `mash sketch -i`, writes `genome_manifest.json` and
`genome_qc.tsv`, and registers the resource only after validation. The full
database is intentionally not bundled with this repository.

## What is computed

Mash `dist` ranks nearest reference phages by Mash distance, p-value, and
matching hashes. Mash is a screening stage; it is not ANI and is not formal
taxonomy. For zero-distance hits, PhageMine may confirm identical,
reverse-complement, or rotation-equivalent nucleotide sequence. If `blastn`
and reference sequences are available, PhageMine additionally calculates its
own bidirectional, length-normalized nucleotide similarity and alignment
coverage.

Reference metadata may include accession, description, host genus, phage genus,
subfamily, and family. These fields describe the reference record and do not
assign the query to an ICTV taxon. A 95% species working boundary is reported
as requiring ICTV review; no generic genus boundary is configured, and
partial/poorly aligned queries are withheld from numerical taxonomy.

## Output files

Outputs are written under `comparative/` for a single-genome run (or the
discovery output directory for batch workflows):

- `inphared_nearest_phages.tsv` and `.json` retain ranked accession-level hits;
- `inphared_summary.tsv` collapses duplicate reference accessions representing
  the same described phage;
- `run_manifest.json` records resource version, manifest path, commands, and
  interpretation; `report.html`/`report.md` include a nearest-reference table.

Important columns are `reference_accession`, `mash_distance`, `p_value`,
`matching_hashes`, `reference_description`, host and phage metadata,
`intergenomic_similarity_percent`, query/reference alignment percentages,
`genome_length_ratio`, `similarity_method`, taxonomy eligibility and scope,
sequence confirmation, and interpretation.

If the resource is not registered, fails validation, or Mash is unavailable,
the output records `status=SKIPPED` with a reason. This is distinct from a
completed comparison with no matching references. A missing `blastn` leaves
Mash screening available while similarity and taxonomic interpretation remain
not calculated.

## Limitations

INPHARED results are comparative context, not proof of common ancestry,
species identity, host range, or ICTV classification. Similarity requires
appropriate alignment coverage and independent taxonomic review. INPHARED
product labels and metadata are computational/reference annotations. The
PhageMine v1.2 seven-phage structural benchmark did not validate INPHARED
performance; that benchmark should not be presented as an INPHARED accuracy
study.
