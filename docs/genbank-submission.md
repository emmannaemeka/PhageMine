# GenBank pre-submission workflow

PhageMine produces a **local pre-submission package**. PhageMine QC, PhageMine pre-submission validation, optional NCBI `table2asn` validation, and final NCBI submission/review are separate stages. No PhageMine command contacts NCBI.

Provide real user metadata with `--metadata metadata.json`. Missing required fields produce `INCOMPLETE_METADATA`; PhageMine never invents them. Required fields are organism, phage name, host, isolation source, collection date, geographic location, sequencing technology, submitter name/email, and authors.

If `table2asn` is on `PATH`, PhageMine runs it locally and preserves `.val`, `.stats`, and other generated outputs in `genbank_submission/table2asn_output`. Otherwise `validation.json` explicitly says official NCBI validation was not performed and includes installation guidance. Obtain the official `table2asn` distribution and current instructions from NCBI's Submission Portal tools documentation, install it locally, and add it to `PATH` (or pass `--table2asn /path/to/table2asn`). Review all output with the current NCBI requirements before submission.

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
