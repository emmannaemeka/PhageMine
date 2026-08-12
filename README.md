# PhageMine

PhageMine is an evidence-first prototype for annotating bacteriophage genomes and prioritising uncharacterised proteins for experimental follow-up.

> Computational evidence produces hypotheses; it does not establish biological function.

## Quick start

The dependency-light demonstration uses a built-in ORF caller and a deliberately labelled mock evidence backend:

```bash
PYTHONPATH=src python -m phagemine run examples/demo_phage.fasta --output results/demo
```

It writes GFF3, CDS/protein FASTA, annotation TSV, evidence JSON, a candidate ranking TSV, Markdown/HTML reports, and a reproducibility manifest. Run `python -m unittest discover -s tests` for tests.

`phagemine genbank genome.fasta` also runs the local-only pre-submission package builder. It creates `genome.fsa`, `features.tbl`, `proteins.faa`, `validation.json`, provenance, and review instructions under `genbank_submission/`. It does not contact NCBI or submit anything.

Pass real submission information with `--metadata metadata.json`; PhageMine never fills in absent source or submitter fields. The future GUI should collect the same metadata before it enables final package generation. See `docs/genbank-submission.md` for the separation between PhageMine checks, optional `table2asn` validation, and NCBI's final review.

## Scientific scope

The three scores are deliberately independent:

- **Functional confidence**: strength of support for an interpretation.
- **Biological interest**: priority for investigation, not probability of function.
- **Evidence diversity**: number of distinct evidence modalities supporting priority.

The built-in mock backend is only for exercising the complete workflow without external databases. It is marked in every output and must be replaced by real similarity/domain/conservation adapters for biological use. See `docs/evidence-framework.md` and `docs/database-policy.md`.
