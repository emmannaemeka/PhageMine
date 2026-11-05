# GenBank pre-submission workflow

PhageMine produces a **local pre-submission package**. PhageMine QC, PhageMine pre-submission validation, optional NCBI `table2asn` validation, and final NCBI submission/review are separate stages. No PhageMine command contacts NCBI.

Provide real user metadata with `--metadata metadata.json`. Missing required fields produce `INCOMPLETE_METADATA`; PhageMine never invents them. Required fields are organism, phage name, host, isolation source, collection date, geographic location, sequencing technology, submitter name/email, authors, and the NCBI-assigned locus tag prefix.

If `table2asn` is on `PATH`, PhageMine runs it locally and preserves `.val`, `.stats`, and other generated outputs in a fresh `genbank_submission/table2asn_*` directory. Otherwise `validation.json` explicitly says official NCBI validation was not performed and includes installation guidance. Obtain the official `table2asn` distribution and current instructions from NCBI's Submission Portal tools documentation, install it locally, and add it to `PATH` (or pass `--table2asn /path/to/table2asn`). Review all output with the current NCBI requirements before submission.

## Curator corrections and revision

Metadata and annotation corrections can be applied without rerunning gene
prediction or evidence searches:

```bash
phagemine revise results/run --corrections submission_corrections.yaml \
  --output results/run/revised_submission
```

The original computational annotation remains unchanged. The revised package
records original and corrected values, reason, source, timestamp, sequence
checksum, parent submission version, and correction-file checksum. A changed
nucleotide FASTA requires a new PhageMine analysis because gene calls and
downstream evidence may change. Internal validation is local pre-submission
validation only; it is not NCBI acceptance.

## Annotation export and validation status

Each CDS is accompanied by a gene feature with the same locus tag. Supported final
computational product annotations are preserved; weak, hypothetical and mock claims
remain hypothetical proteins. The package includes both `features.tbl` for inspection
and `genome.tbl`, whose basename matches `genome.fsa` for table2asn. RNA additions
update both tables and rerun official validation.

Official validation checks `.val` findings and fatal `.dr` discrepancies, not just the
process exit code. The report distinguishes `READY`, `READY_WITH_WARNINGS`,
`INCOMPLETE_METADATA`, `REQUIRES_REVIEW`, `VALIDATION_FAILED`, and
`PENDING_NCBI_VALIDATION`. No package is marked ready when official validation is
unavailable or fails. Local readiness does not establish NCBI acceptance.

Use an NCBI-generated submission template containing your actual author and
affiliation information:

```bash
export PHAGEMINE_SUBMISSION_TEMPLATE=/absolute/path/to/template.sbt
```

The generated minimal template is syntactically readable, but missing affiliation
fields can still fail official validation. Supply source metadata and the approved
`locus_tag_prefix`; PhageMine does not invent missing values or alter genome bases
to clear errors. Biological boundary decisions, completeness, and genetic code
require review against evidence. Consult [NCBI validation explanations](https://www.ncbi.nlm.nih.gov/genbank/genome_validation/)
and [feature-table formatting](https://www.ncbi.nlm.nih.gov/genbank/feature_table/).
