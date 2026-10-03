# Researcher readiness and adoption criteria

Assessment date: 2026-10-03. Software execution and scientific accuracy are
separate questions. PhageMine remains beta research software. It supports
exploratory annotation and review; broad superiority has not been demonstrated.

## Why a researcher might decline it

| Concern | Evidence or current behavior | Resolution |
|---|---|---|
| Functional accuracy is uncertain | Automated adjudication left 476/689 units unresolved; named-product yield is not precision | Independent blinded expert review on held-out genomes remains required |
| Structural agreement is weaker on the small panel | Strict F1 0.7762 versus Prokka 0.8846 on seven curated genomes | Review disagreements; do not optimize against this test panel |
| A first run is difficult | Several external tools and large downloads | Core-first onboarding, explicit resource choices, PHANOTATE compatibility dependency |
| Installation might succeed without execution | Historical release check skipped example execution | Built-wheel execution outside checkout and real-caller Core smoke checks |
| The GUI can misrepresent results | Early manifests treated as complete, Both results missed, protein IDs ambiguous | Conservative manifest checks, native layout resolution, sample-qualified selection |
| Complete genome annotation is needed | Structured RNAs and exceptional gene features are not automatically called | Specialist review remains required; see limitations |
| Operational scale is unclear | No validated broad runtime/memory comparison | Resource-backed scale tests and comparable timing provenance remain required |

## Changes in this working branch

The GUI defaults to localhost, validates input types and duplicate upload names,
expands home paths, avoids duplicate Both annotations, loads discovery results
from their actual directory, distinguishes samples in protein selection,
reports malformed artifacts and excludes symlinks from exports. Submission
FASTA is included in downloads. ZIP construction happens only on request.
Operational Doctor checks and provider readiness are visible in the GUI.

Batch commands return a failing exit status if any sample fails. Empty batches
are rejected, and Both stops before discovery when annotation has failures,
preserving completed outputs and failure records for recovery.

CI installs a built wheel into a separate environment, checks dependencies and
executes it outside the repository with a disposable registry. GUI tests cover
page rendering and successful runs. External-tool CI runs a real-caller Core
smoke check on relevant pull requests and scheduled/manual runs. The source
distribution includes the environment, examples and smoke-check script.

The environment retains `setuptools<81` for PHANOTATE 1.6.7's `pkg_resources`
import. Full evidence installation and platform-specific Conda availability
must still be checked separately. See [current software validation](../RESEARCHER_VALIDATION.md)
for actual local results. Adding a workflow is not a claim that remote CI passed.

## Evidence required before broad recommendation

Follow the [validation protocol](validation-protocol.md) and
[independent evidence-engine protocol](../validation/independent_evidence_engine/README.md).
Freeze development/test partitions, references and database snapshots. Identify
database overlap and report a leakage-controlled held-out subset. Independent
experts must adjudicate unresolved units without seeing tool labels and retain
disagreement and abstention records.

Report coordinate precision/recall, reviewed product assertion precision,
coverage, unsupported specificity and abstention with uncertainty at the genome
level. Assess each evidence-strength category separately; rule-based labels
cannot be called calibrated probabilities without results. Include diverse
genomes and separate challenge sets. Preserve the historical benchmark and
publish new results as a distinct evaluation.

Before release, run full resource installation and deep checks, full-evidence
annotation, a discovery/Both cohort and resume recovery using the published
artifact. Record tool/database versions, output integrity, runtime and peak
memory. Exercise Linux and supported macOS architectures. Have external
researchers attempt installation, inspect uncertain records and export results;
record task completion, mistakes and feedback.

These activities require reference data, external tools, computational resources
and independent reviewers. Code changes cannot substitute for those results.
Retain the beta designation until the relevant [release gates](release-checklist.md)
are satisfied.
