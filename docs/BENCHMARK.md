# PhageMine v1.2 benchmark

The corrected benchmark evaluates seven curated reference phages (T4, Lambda,
T7, T5, PhiX174, P22, and Mu). The full local evidence is preserved separately;
this repository distributes only compact summaries, provenance, figures, and
checksums under [`docs/benchmark_v1.2/`](benchmark_v1.2/).

## Measured results

| Tool | CDS | Strict F1 | Relaxed F1 | Named products |
|---|---:|---:|---:|---:|
| PhageMine | 804 | 0.7762 | 0.8821 | 439 |
| Pharokka | 804 | 0.7762 | 0.8821 | 398 |
| Prokka | 675 | 0.8846 | 0.9352 | 328 |

PhageMine and Pharokka had identical CDS coordinates across all seven genomes
in this configuration. Their structural scores are therefore identical; their
functional outputs differ downstream. Prokka showed stronger reference-relative
structural agreement on this curated panel.

Strict scoring requires exact accession, start, stop, and strand. Relaxed scoring
counts a same-strand prediction as detecting a reference CDS when it covers at
least 50% of the reference span and has at least 20% reciprocal overlap. Of 216
PhageMine strict non-exact calls, 141 were alternative-boundary or same-strand
overlap cases. These should not automatically be interpreted as false
biological genes. Non-reference calls remain non-reference predictions, not
proven novel CDSs.

![Strict versus relaxed F1](benchmark_v1.2/figures/strict_vs_relaxed_f1.png)

![Per-genome structural performance](benchmark_v1.2/figures/per_genome_strict_f1.png)

![Functional annotation yield](benchmark_v1.2/figures/functional_yield.png)

## Reproducing the summary

The source tables are in `docs/benchmark_v1.2/tables/`. The included figures
were generated from those tables; no biological databases are required to view
them. To regenerate equivalent plots, run:

```bash
python docs/benchmark_v1.2/scripts/make_figures.py
```

Named-product yield is a lexical output metric, not functional accuracy. The
seven-genome panel is limited, exact coordinates penalize alternative boundaries,
and runtime ranking is unsupported. The historical T4 PHANOTATE 294-versus-297
Linux discrepancy remains unresolved. No universal tool ranking is claimed.
