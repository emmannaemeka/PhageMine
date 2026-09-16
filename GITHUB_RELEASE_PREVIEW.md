# GitHub release preview — PhageMine v1.2

This is the reviewed release candidate for publication. Commit and push actions
are being performed only after the final audit in this document.

## Repository audit

- Branch: `v1.2-finalization`
- Tracking branch: `github/v1.2-finalization`
- Remotes: `github` points to the public GitHub URL; `origin` is a local
  publication-workspace path and is not release metadata.
- Initial working tree size: approximately **3.8 GiB** across 12,694 files,
  with 140 files tracked before this candidate.
- Files above 10 MiB, 50 MiB, and 100 MiB (50, 10, and 4 files respectively)
  are generated benchmark JSON/FASTA
  or intermediate outputs under local evaluation trees. None are intended for
  GitHub.
- No credentials, API keys, passwords, or private tokens were found in
  tracked source/documentation. Personal absolute paths occur in local,
  ignored development captures and are excluded from the candidate.
- `.DS_Store`, Python caches, pytest caches, logs, temporary directories,
  local environments, databases, MMseqs/DIAMOND indexes, and raw benchmark
  evidence are ignored. They are not deleted.

## Files to add or modify

The candidate changes release metadata and documentation (`README.md`,
`CHANGELOG.md`, `CITATION.cff`, `LICENSE` retained,
`pyproject.toml`, `docs/installation.md`, `docs/INPHARED.md`,
`docs/BENCHMARK.md`, `INPHARED_FEATURE_AUDIT.md`, and this validation/preview
documentation), adds the small-genome Prodigal mode regression test and the
validated-output export module, and carries the existing reporting integration.
The public compact benchmark package is under `docs/benchmark_v1.2/` and
contains selected TSV tables, PNG/SVG/PDF figures, provenance, and a portable
figure regeneration script.

## Intentionally excluded

The full `evaluation/v1.2_tool_comparison/` tree, frozen original and corrected
benchmark evidence, final structural interpretation workspace, raw outputs,
logs, status captures, downloaded databases, local Conda environments,
diagnostic runs, temporary files, and machine-specific development transcripts
remain on disk but are excluded from GitHub. Third-party INPHARED, PHROGs,
Pfam, VOGDB, Swiss-Prot, PMFDB, and INPHARED-derived databases are acquired by
the documented installers and are not redistributed.

## Benchmark material included

`docs/BENCHMARK.md` reports the frozen corrected seven-phage structural and
functional summary. It states that PhageMine and Pharokka have identical CDS
coordinates in this panel, that Prokka has stronger reference-relative
structural metrics, and that named-product yield is not functional accuracy.
The compact package includes strict/relaxed metrics, category counts, reference
provenance, count differences, functional yields, figures, and checksums.

## INPHARED status

The feature audit classifies INPHARED as **A: fully implemented and currently
callable**. It is automatically consumed by `annotate`, `run`, and batch
`discover`/`both` workflows after validation. It uses Mash screening followed
by optional bidirectional BLASTN length-normalized similarity, reports reference
metadata, and does not make formal ICTV assignments. See `docs/INPHARED.md` and
`INPHARED_FEATURE_AUDIT.md`.

## Validation status

The clean installed clone imported version 1.2.0, displayed CLI help, and
passed 272 tests with 10 deselected. The INPHARED/resource/small-genome subset
passed 31 tests. Details, including the offline pip build-isolation note, are
in `RELEASE_VALIDATION.md`. No annotation benchmark was rerun.

## Review decisions recorded

The reviewed candidate retains the existing author and MIT license metadata;
it invents no DOI, ORCID, affiliation, or new license terms. The additive
validation export and small-genome regression are included, and compact
benchmark summaries are distributed while full evidence remains local. The
existing `github` remote and `v1.2-finalization` branch are the publication
targets.

## Release gate

The candidate is ready for commit, push, tag, and release publication.
