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

## When NCBI returns a submission with errors

A submitted genome can be edited only while its status in the GenBank
Submission Portal is **Error**. Press **Fix**; if there is no Fix button,
click the submission's SUB# and use the Sequencing Technology, Sequences,
Source Modifiers or References tabs. Paste NCBI's messages into Ask PhageMine:
it sorts them into fixable wording, missing source details, and problems
needing gene-model review, and prepares a revised package after you approve.
Upload the **complete** revised package, every sequence, never only the
corrected one, then click **Submit**; the status stays Error until you do.
NCBI explains individual errors at
[genbank/validation](https://www.ncbi.nlm.nih.gov/genbank/validation/) and,
for viruses, [genbank/sequencecheck/virus](https://www.ncbi.nlm.nih.gov/genbank/sequencecheck/virus/).

## NCBI requirements review

Every package is reviewed against published NCBI and INSDC requirements
(`phagemine/submission_review.py`). Each finding in `validation.json` and
`NCBI_VALIDATION_REPORT.txt` names the NCBI check it mirrors, the fix, and
the NCBI page it comes from. ERROR findings are expected to be flagged or
returned by NCBI and make the package `VALIDATION_FAILED`; WARNING findings
are review items. `table2asn` and GenBank staff remain the authority.

| Area | What is checked | NCBI check |
| --- | --- | --- |
| Collection date | DD-Mon-YYYY, Mon-YYYY, YYYY or ISO 8601; real and not in the future; or an INSDC missing-value term. Required on all submissions since December 2024 | SEQ_DESCR_BadCollectionDate |
| Location | Country/ocean on the INSDC geo_loc_name list, then a colon and the region (`Nigeria: Plateau State, Jos`). Required since December 2024 | SEQ_DESCR_BadGeoLocNameCode |
| Organism | `<host genus> phage <name>` form (warning) | taxonomy review |
| Sequence | Under 200 bp, terminal Ns, runs of 100+ Ns, invalid characters, high ambiguity, a missing base | SHORT_SEQUENCES, TerminalNs, InternalNsInSeqRaw, DISC_PERCENT_N, ZERO_BASECOUNT |
| Product names | Exactly "hypothetical protein" for unknown function; "putative", not predicted/probable/possible; no homolog, database IDs, locus tags, EC numbers, partial/fragment, trailing periods or forbidden punctuation; lowercase, American spelling | SUSPECT_PRODUCT_NAMES |
| CDS | Partial only at a sequence end; complete CDSs start with a translation-table-11 start codon; proteins under 50 aa; contained CDSs; overlapping CDSs with the same name | PARTIAL_CDS_COMPLETE_SEQUENCE, SEQ_FEAT.StartCodon, SHORT_PROT_SEQUENCES, CONTAINED_CDS, OVERLAPPING_CDS |
| Locus tags | Prefix 3-12 letters/digits starting with a letter; written on gene features only | BAD_LOCUS_TAG_FORMAT, MISSING_GENES |
| Template | PhageMine's template has names only; supply NCBI's official template | GENERIC.MissingPubInfo |

What the package now contains for NCBI:

- `genome.fsa` carries the source information on the definition line, where
  `table2asn` reads it: organism, host, isolation source, collection date,
  geo_loc_name, isolate, strain, molecule type, and `[gcode=11]` (bacteria and
  their phages).
- `features.tbl` has a `gene` feature (with the locus tag) and a `CDS` for
  every protein, and a `gene` for every tRNA. No local file references are
  written into notes.
- `proteins.faa` gives complete CDSs that start with an alternative start
  codon (GTG, TTG...) an initial methionine, as NCBI translates them.

Recommended metadata additions:

```json
{
  "organism": "Xanthomonas phage vB_XcaM_JOS1",
  "collection_date": "14-Mar-2025",
  "geographic_location": "Nigeria: Plateau State, Jos",
  "locus_tag_prefix": "JOS1",
  "submission_template": "template.sbt"
}
```

Register a BioProject to obtain the locus tag prefix, and create
`template.sbt` with NCBI's Submission Template form.

Sources: [Common discrepancy reports](https://www.ncbi.nlm.nih.gov/genbank/asndisc.examples/),
[Validation error explanations](https://www.ncbi.nlm.nih.gov/genbank/genome_validation/),
[International protein nomenclature guidelines](https://www.ncbi.nlm.nih.gov/genbank/internatprot_nomenguide/),
[Genome annotation guide](https://www.ncbi.nlm.nih.gov/genbank/genomesubmit_annotation/),
[Source modifiers](https://www.ncbi.nlm.nih.gov/WebSub/html/help/genbank-source-table.html),
[INSDC geo_loc_name list](https://www.ncbi.nlm.nih.gov/genbank/collab/country/),
[INSDC missing values](https://www.insdc.org/technical-specifications/missing-value-reporting/),
[Sequence checks](https://www.ncbi.nlm.nih.gov/genbank/sequencecheck/virus/).
