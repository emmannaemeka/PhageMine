# Changelog

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
