# Changelog

## 1.0.6

- Separate domain evidence notes from final product names; organism-specific
  Pfam descriptions remain in raw evidence and cannot become phage products.
- Label functional strength as rule-based and explicitly uncalibrated.
- Remove the universal 70% genus rule, abstain from numerical-taxonomy
  interpretation for partial/poorly aligned sequences, and stop describing the
  internal BLASTN calculation as VIRIDIC output.
- Mark hallmark results as annotation-text screening rather than independent
  profile validation or genome-completeness estimates.
- Add truth-set benchmark inputs and calibration-ready performance outputs.
- Add reproducible environment metadata, declared plotting/test dependencies,
  package metadata, and automated multi-platform checks.

## 1.0.5

- Isolate tests from the persistent user database registry and make registry
  writes atomic with a recoverable backup.
- Add `phagemine doctor --deep` to open and validate the registered PHROGs
  MMseqs2 and PyHMMER databases, rather than inferring readiness from filenames.
- Report MMseqs2 and PyHMMER PHROGs backend status in normal run progress.
- Promote repeated, concordant, high-coverage annotated PHROGs profile hits
  over a single unannotated profile and add Ijeoma-derived HNH/terminase
  regression cases.
- Keep functional uncertainty separate from gene-call confidence; exact
  PHANOTATE/Prodigal coordinate and strand agreement is now reported
  consistently and is never converted into a possible-false-ORF flag.
- Add gene-call confidence and independent functional/gene-call review fields
  to researcher-facing annotation tables.
- Propagate `--threads` from `phagemine batch` into every single-genome run and
  repair nested progress-stage reporting.

## 1.0.4

- Add sensitive PHROGs profile-HMM searching through PyHMMER while retaining
  the existing MMseqs2 search.
- Extract and register the existing `all_phrogs.h3m` from the checksum-pinned
  Pharokka PHROGs distribution; no additional database download is required.
- Deduplicate identical PHROG calls across MMseqs2 and PyHMMER so one
  biological source cannot be counted twice, while recording backend
  corroboration and search-specific statistics.
- Add `phagemine databases attach-phrogs-hmm` for attaching an already
  installed Pharokka `all_phrogs.h3m` to an existing PHROGs registration.
- Extend PHROGs refresh/resume to replace older PHROGs evidence rather than
  silently duplicating it, and allow targeted profile-HMM refreshes without
  rerunning gene prediction or unrelated database searches.
- Report PyHMMER and sensitive PHROGs profile-search readiness explicitly in
  `phagemine doctor`.

## 1.0.3

- Replace the invalid `1 - Mash distance` pseudo-similarity with Mash-only
  candidate screening followed by VIRIDIC-compatible bidirectional BLASTN
  whole-genome intergenomic similarity.
- Report query/reference aligned fractions, genome-length ratio, explicit
  species/genus thresholds, threshold provenance, and cautious numerical-
  taxonomy interpretations.
- Add the completed INPHARED numerical-taxonomy table and downloads to the
  primary `report.html` after comparative discovery finishes.
- Add simple Prokka-style `gene`, `product`, `EC_number`, evidence-source and
  biotechnology-relevance columns. Gene symbols and EC numbers are transferred
  only from explicit, strong, non-conflicting evidence.
- Normalize defensible phage terms including major capsid, Hoc-like head
  decoration, head-scaffolding, terminase, endolysin and helicase products;
  suppress unsupported fungal, apoptosis and bacterial-envelope labels.
- Identify Hoc-like decoration proteins as potential capsid-display candidates
  in a separate application note that explicitly requires experimental
  confirmation.
- Correct the progress denominator for runs that include PHANOTATE/Prodigal
  reconciliation and alternative-ORF adjudication.
- Allow Core `run` to complete without invoking MMseqs2 discovery when no
  comparative resources are registered.

## 1.0.2

- Simplify the primary annotation table to protein ID, coordinates, strand,
  length, plain-language classification, proposed function, confidence, best
  evidence, and a review flag. Detailed evidence remains in the dedicated
  classification and protein-detail outputs.
- Add `annotated_proteins.faa` with the proposed product, coordinates, strand,
  and confidence in every FASTA header.
- Make full single-genome `run` perform PHANOTATE/Prodigal comparison and write
  explicit gene-call confidence and review tables without silently deleting
  caller-specific ORFs.
- Place results in `./<input stem>_phagemine_results` when `--output` is not
  supplied; an explicit output path continues to take precedence.
- Add a duplicate-collapsed INPHARED summary and distinguish exact sketch
  matches from nearest-neighbour screens while retaining the accession-level
  audit table.
- Extend annotation comparison imports to Pharokka, Phold, multiPhATE2, and
  Prokka, with separate gene-model and product-name comparisons.
- Add an explicit six-tier evidence hierarchy to the detailed classification
  output; the hierarchy contains no structural-search dependency.
- Add conservative hallmark-system checks and a single prioritized
  `annotation_review.tsv` combining function conflicts, disputed gene calls,
  and hallmark components not established by current evidence.
- Confirm zero-distance INPHARED matches by built-in exact nucleotide,
  reverse-complement, and rotation-equivalence checks when the reference
  sequence is available; do not infer ANI or taxonomy from Mash.
- Record checksummed external-tool inputs and hallmark-name comparisons in the
  reproducible annotation benchmark output.
- Do not infer or report phage lifestyle.

## 1.0.1

- Synchronize the release version with the post-1.0.0 database installer and
  INPHARED functionality published on GitHub.
- Make executable validation fail closed: a binary is `READY` only when its
  version/help probe exits successfully. Dynamic-linker failures and other
  runtime errors are now reported as `BROKEN` with their diagnostic text.
- Make `phagemine doctor` return a non-zero process status when a required
  production executable is missing or broken; Prodigal and table2asn remain
  optional for the normal annotation and pre-submission workflows.
- Distinguish an installed-but-blocked database from a database that is not
  installed, and recommend executable repair instead of unnecessary database
  reinstallation.
- Load unknown future registry resource types without crashing access to all
  existing registered databases.
- Require Mash in the Bioconda run dependencies for whole-genome INPHARED
  comparison.

- Add `phagemine databases install pfam`, `vogdb`, `swissprot`, `phrogs`,
  `pmfdb`, and `inphared`.
- Add `phagemine databases install --all` with resumable downloads, published
  checksum verification, database preparation, manifests, registration, and
  fail-closed readiness validation.
- Add post-install database guidance to top-level help, the Bioconda post-link
  template, and `phagemine doctor` recommendations.
- Build PMFDB and the whole-genome Mash resource from the checksum-pinned
  INPHARED 7 April 2026 release, with conservative annotation semantics,
  resource QC, provenance, automatic discovery use, and nearest-phage outputs.

## 1.0.0

The first stable public release. This release provides evidence-based
single-genome annotation, batch annotation, Discovery Mode, and Both Mode;
PHANOTATE integration; Pfam/VOGDB/Swiss-Prot/PHROGs evidence; conservative
evidence fusion and classification; genomic context/modules; candidate
ranking; PMF clustering and recurrence; PMFDB validation; discovery ranking;
publication figures with source data; HTML reports; protein/CDS retrieval;
GenBank pre-submission packaging; resume/checkpointing; and evidence reuse
with provenance.

The release does not bundle external evidence databases, PMFDB, or research
datasets. Unknown and no-match states are not interpreted as novelty.

## 1.0.0rc1 (historical)

This release candidate freezes the v1.0 computational workflow. It provides
the lightweight PhageMine Core FASTA-to-report workflow, optional evidence
integration, PMF/PMFDB integration when explicitly requested, local GenBank
pre-submission validation, and submission revision from an explicit correction
file without rerunning biological analysis.

The release does not bundle external evidence databases, PMFDB, or research
datasets. See `docs/limitations.md` for interpretation and validation limits.
