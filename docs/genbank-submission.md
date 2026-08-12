# GenBank pre-submission workflow

PhageMine produces a **local pre-submission package**. PhageMine QC, PhageMine pre-submission validation, optional NCBI `table2asn` validation, and final NCBI submission/review are separate stages. No PhageMine command contacts NCBI.

Provide real user metadata with `--metadata metadata.json`. Missing required fields produce `INCOMPLETE_METADATA`; PhageMine never invents them. Required fields are organism, phage name, host, isolation source, collection date, geographic location, sequencing technology, submitter name/email, and authors.

If `table2asn` is on `PATH`, PhageMine runs it locally and preserves `.val`, `.stats`, and other generated outputs in `genbank_submission/table2asn_output`. Otherwise `validation.json` explicitly says official NCBI validation was not performed and includes installation guidance. Obtain the official `table2asn` distribution and current instructions from NCBI's Submission Portal tools documentation, install it locally, and add it to `PATH` (or pass `--table2asn /path/to/table2asn`). Review all output with the current NCBI requirements before submission.
