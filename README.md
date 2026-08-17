# PhageMine

PhageMine is an evidence-first prototype for annotating bacteriophage genomes and prioritising uncharacterised proteins for experimental follow-up.

> Computational evidence produces hypotheses; it does not establish biological function.

## Quick start

### Installation

Build and install the release candidate with `python -m build` followed by
`python -m pip install dist/phagemine-1.0.0rc1-py3-none-any.whl`. Production
Core runs require Python 3.10+ and an installed PHANOTATE executable/script.

```bash
phagemine run genome.fasta --phanotate /path/to/phanotate.py --output results/run
```

Core performs genome QC, PHANOTATE gene prediction, conservative reporting,
candidate prioritization, and local GenBank pre-submission packaging. Optional
Pfam/HMMER, VOGDB/HMMER, PHROGs/MMseqs2, Swiss-Prot/DIAMOND, and Prodigal
resources are enabled explicitly or through the local resource registry; they
are not required for a Core installation and are never downloaded automatically.
PMF family/PMFDB workflows are separate optional workflows. GenBank output is a
local pre-submission package; PhageMine does not claim NCBI acceptance. Results
are evidence-based hypotheses and the software makes no automatic novelty claim.

The dependency-light demonstration uses a built-in ORF caller and a deliberately labelled mock evidence backend:

```bash
# Production default: locally installed PHANOTATE
PYTHONPATH=src python -m phagemine run genome.fasta --output results/run

# Synthetic fixture only (no PHANOTATE required)
PYTHONPATH=src python -m phagemine run examples/demo_phage.fasta --gene-predictor demo --output results/demo
```

It writes GFF3, CDS/protein FASTA, annotation TSV, evidence JSON, a candidate ranking TSV, Markdown/HTML reports, and a reproducibility manifest. Run `python -m unittest discover -s tests` for tests.

### Test suites

Routine development tests exclude resource-backed integration tests:

```bash
PYTHONPATH=src python -m pytest -q -m "not integration"
```

The resource-backed tests remain available explicitly:

```bash
PYTHONPATH=src python -m pytest -q -m integration
```

The complete suite (overriding the default marker filter) is
`PYTHONPATH=src python -m pytest -q -o addopts=''`.
Integration tests retain coverage for registered HMMER, DIAMOND, MMseqs2, and
end-to-end resource-backed runs; they are excluded from the default command
only to keep routine regression feedback fast.

`phagemine genbank genome.fasta` also runs the local-only pre-submission package builder. It creates `genome.fsa`, `features.tbl`, `proteins.faa`, `validation.json`, provenance, and review instructions under `genbank_submission/`. It does not contact NCBI or submit anything.

Pass real submission information with `--metadata metadata.json`; PhageMine never fills in absent source or submitter fields. The future GUI should collect the same metadata before it enables final package generation. See `docs/genbank-submission.md` for the separation between PhageMine checks, optional `table2asn` validation, and NCBI's final review.

### Revising a submission package

After curator or researcher review, metadata and explicit feature corrections can
be applied without rerunning biological analysis:

```bash
PYTHONPATH=src python -m phagemine revise RESULTS_DIR \
  --corrections submission_corrections.json \
  --output revised_submission
```

The correction file is a transparent JSON object with `metadata` and
`annotations` sections. Annotation entries may set `product`, `note`, `partial`,
or `locus_tag`, and should include `reason` and `source`. The original analysis
files are never overwritten; the revised package is written under
`submission_v2/` with correction audit TSV/JSON and sequence/checksum
provenance. A changed genome checksum fails with “Genome sequence has changed; a
new PhageMine analysis is required.” Metadata/annotation corrections revise only
the submission package; sequence corrections require a new PhageMine run.

## Scientific scope

The three scores are deliberately independent:

- **Functional confidence**: strength of support for an interpretation.
- **Biological interest**: priority for investigation, not probability of function.
- **Evidence diversity**: number of distinct evidence modalities supporting priority.

The built-in mock backend is only for exercising the complete workflow without external databases. It is marked in every output and must be replaced by real similarity/domain/conservation adapters for biological use. See `docs/evidence-framework.md` and `docs/database-policy.md`.

## Gene prediction

Core production analysis uses the `PHANOTATEPredictor` adapter by default. Install PHANOTATE separately and ensure `phanotate.py` or `phanotate` is on `PATH`, or use `--phanotate /path/to/phanotate.py`. It is intentionally not bundled or installed by PhageMine. The predictor identity, executable, version, and parameters are recorded in the run manifest. The legacy simple ORF caller is retained only as `--gene-predictor demo` for test fixtures.

## Sequencing provenance

PhageMine analyzes assembled FASTA independently of sequencing platform. Optionally pass `--sequencing-provenance provenance.json` to record platform, library, assembly, polishing, and raw-read availability. Without it, the platform is explicitly `UNKNOWN`; PhageMine never infers sequencing technology from FASTA. Sequencing provenance is separate from genome topology, orientation, and rotation.

## Optional Pfam evidence

PhageMine can use a user-supplied local Pfam HMM database with `--pfam /path/to/Pfam-A.hmm`. HMMER must be installed separately (or supplied with `--pfam-hmmscan`). PhageMine does not download Pfam. Missing Pfam/HMMER is reported as `UNAVAILABLE` and produces no fabricated domain evidence. Pfam hits are evidence records, not automatic functional assignments.
