# Limitations and interpretation

PhageMine is beta research software. Its outputs are computational,
evidence-supported hypotheses and must not be represented as experimentally
validated functions, formal ICTV assignments, proof of genome completeness, or
automatic NCBI acceptance.

## Functional annotation

Pfam, VOGDB, PHROGs and Swiss-Prot matches can support homology or domain
architecture, but do not by themselves prove biochemical function. Domain
descriptions are retained as evidence notes and are not final product names.
The `HIGH`, `MODERATE` and `LOW` labels are deterministic rule-based evidence
strength categories. They have not yet been calibrated as probabilities on a
large expert-reviewed truth set.

PMFDB `NO_EXTERNAL_MATCH` is not proof of novelty. Searches capped with
`--max-seqs 100` are not prevalence estimates. Database product labels remain
computational predictions unless explicit experimental provenance is present.

## Gene models

PHANOTATE is the retained primary caller. Prodigal comparison measures caller
agreement and creates a review queue; it does not establish which discordant
model is biologically correct. PhageMine does not currently call tRNA, tmRNA,
other structured RNAs, programmed frameshifts, translational bypasses, introns,
or alternative genetic codes. Submitters must review these features with
specialist tools.

## Hallmarks and completeness

`hallmark_completeness.tsv` is retained for backward compatibility, but its
current method is an annotation-text screen. It is not an independent profile-
HMM search and is not a genome-completeness estimator. `NOT_ESTABLISHED` is not
evidence of biological absence.

## INPHARED and taxonomy

Mash selects candidate references only. PhageMine's reported nucleotide value
is its own bidirectional BLASTN length-normalized calculation; it is not
presented as VIRIDIC output. PhageMine does not assign taxa. It withholds
boundary interpretation for partial or poorly aligned queries, applies no
universal genus threshold, and requires users to consult current family-
specific ICTV criteria and formal phylogenetic analyses.

## Benchmarking

Agreement with Pharokka, Phold, multiPhATE2 or Prokka is not truth, especially
when tools share PHANOTATE, PHROGs or other databases. Accuracy or superiority
claims require an expert-reviewed truth set supplied to `phagemine benchmark`.

## Submission

PhageMine GenBank checks are pre-submission checks and do not claim NCBI
acceptance. A nucleotide sequence change requires a new analysis. Optional
`table2asn` validation and final human review remain necessary.
