# PhageMine DeepOrigins — evidence audit prototype

**Status:** experimental, not validated for biological discovery.

## Question

Can apparently unrelated proteins exhibit structurally conserved architectures
across phages, bacteria and eukaryotes? Such similarities can arise from distant
common ancestry, convergence, or misleading comparisons. DeepOrigins must not
treat structural similarity as proof of ancestry or function.

## Implemented in this branch

A small, dependency-free Python evidence auditor accepts normalized TSV observations
from external sequence/structure comparison tools, checks their numeric ranges,
groups evidence by query, and writes an explicit uncertainty-preserving JSON report.
It does **not** run Foldseek, MMseqs2, HMMER, protein folding, or evolutionary
reconstruction. It does **not** change PhageMine's existing annotation pipeline.

Run from an installed checkout:

```bash
python -m phagemine.deeporigins comparisons.tsv --output deeporigins_report.json
```

TSV columns (all required):

```text
query target source_group sequence_identity alignment_coverage structural_similarity method reference
```

Use **tabs** between columns, fractions from 0 to 1 for all three numeric
fields, and source_group = phage, bacteria, eukaryote, or other.
The structural_similarity field is a normalized descriptive value supplied
by the analyst, **not** an interchangeable Foldseek TM-score or probability.
Document how it was computed in the upstream analysis. Reference should
identify a traceable source or dataset version.

## Scientific guardrails

- Review priority is a heuristic queue, not a statistical significance test.
- No thresholds here establish remote homology, ancestry, novelty, or function.
- Compare aligned coverage, fold topology, low-complexity regions, structural
  prediction confidence, reference database redundancy, and negative controls.
- Structural resemblance can reflect convergent evolution.
- Establish independent sequence/profile evidence and genomic-context support.
- Use known homologous and unrelated protein controls before reporting accuracy.
- Maintain provenance, database versions, and held-out benchmark sets.

## Next engineering milestones

1. Map PhageMine protein identifiers to exported protein sequences and annotations.
2. Add adapters for *precomputed* sequence/profile and structural search results.
3. Introduce a reproducible evidence schema with tool/database versions and
   calibrated confidence metrics rather than arbitrary score normalization.
4. Evaluate against curated positive and negative evolutionary-relationship sets.
5. Add a review UI and downloadable provenance-rich candidate dossiers.

No claims of new protein families or ancient ancestry should be made from this
prototype alone.
