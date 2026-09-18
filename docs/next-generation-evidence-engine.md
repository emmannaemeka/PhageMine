# Next-generation evidence engine

The `phagemine.evidence_engine` module implements the evidence-adjudication
layers for development runs. It consumes raw evidence records without
rewriting them, keeps PHANOTATE as the final CDS model, and treats
Pyrodigal-gv as an observational secondary caller. Structural relationships,
structural CDS confidence, family/orthology versus specific-function
confidence, module evidence, architecture hypotheses, architecture-specific
hallmarks, and comparative significance are separate output objects.

`phagemine evidence-engine --input-json PACKET --output DIRECTORY` writes TSV
and JSON provenance artifacts, including `module_evidence`,
`architecture_assessment`, `architecture_hallmarks`, `comparative_assessment`,
`structural_adjudication`, and `functional_adjudication`.

The comparative layer distinguishes a mathematical nearest candidate from a
biologically significant reference. A Mash distance of 1.0 or zero shared
hashes cannot create a genus or species claim. `NO_SIGNIFICANT_REFERENCE` is a
comparative state, not a novelty or taxonomy assignment.

Module status `NOT_ESTABLISHED` is not biological absence. Filamentous-phage
hallmarks are assessed only after the architecture hypothesis; tailed-phage
packaging/tail components are marked `NOT_APPLICABLE_TO_ARCHITECTURE` when that
hypothesis is supported. Replication-initiation evidence is eligible for the
DNA-replication module, while weak subordinate holin/domain evidence cannot
establish lysis against stronger locus-level evidence.

The pinned optional development dependencies are `pyrodigal==3.7.1` and
`pyrodigal-gv==0.3.2`. No biological databases were added. The requested
development FASTA inputs were not present in this checkout, so no genome ID or
expected outcome was fabricated.
