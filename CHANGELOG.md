# Changelog

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
