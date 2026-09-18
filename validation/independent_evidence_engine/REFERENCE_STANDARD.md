# Reference standard

GenBank product text is evidence, not automatic truth. Each reference locus
receives one prespecified reference-standard tier:

- `TIER_1_EXPERIMENTALLY_ESTABLISHED`: direct experiment or primary literature
  establishing the function for the protein or a defensible ortholog.
- `TIER_2_EXPERT_CURATED`: expert-reviewed record with concordant primary or
  reviewed evidence and an informative product.
- `TIER_3_STRONGLY_ORTHOLOGY_SUPPORTED`: full-length, high-coverage,
  biologically coherent orthology with concordant family/context evidence.
- `TIER_4_COMPUTATIONALLY_SUPPORTED`: useful computational evidence without
  the strength required for the primary tiers.
- `UNRESOLVED`: hypothetical, conflicting, or insufficient evidence.
- `NOT_EVALUABLE`: the locus or function cannot be assessed under the frozen
  rules.

The hierarchy is experimental characterization, expert-curated annotation,
reviewed protein records, strongly established orthology, conserved
structural/functional evidence, and literature-supported computational
interpretation. Family membership alone does not establish a specific
biochemical function; a domain hit does not establish the whole-protein
function.

The primary functional analysis uses Tiers 1–3. Tier 4 is included only in a
prespecified expanded-tier sensitivity analysis. At least two independent
reviewers adjudicate difficult semantic cases while blinded to tool identity.
Accepted synonym rules are version-controlled before review. Conflicts, missing
evidence, and reference uncertainty are retained in the reference packet and
excluded or marked unresolved according to the statistical plan.
