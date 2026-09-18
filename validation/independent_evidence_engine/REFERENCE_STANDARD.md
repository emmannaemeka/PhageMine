# Reference standard

GenBank product text is evidence, not automatic truth. Each reference locus
receives a reference-quality state:

- `EXPERIMENTALLY_ESTABLISHED`: direct experimental or literature-supported
  function.
- `STRONGLY_CURATED`: expert-reviewed record with concordant evidence.
- `COMPUTATIONALLY_SUPPORTED`: strong orthology/profile/structure and
  provenance, without direct experimental confirmation.
- `UNRESOLVED`: hypothetical, conflicting, or insufficient evidence.
- `NOT_EVALUABLE`: the locus or function cannot be assessed under the frozen
  rules.

The hierarchy is experimental characterization, expert-curated annotation,
reviewed protein records, strongly established orthology, conserved
structural/functional evidence, and literature-supported computational
interpretation. Family membership alone does not establish a specific
biochemical function; a domain hit does not establish the whole-protein
function.

At least two independent reviewers adjudicate difficult semantic cases while
blinded to tool identity. Accepted synonym rules are version-controlled before
review. Conflicts, missing evidence, and reference uncertainty are retained in
the reference packet and excluded or marked unresolved according to the
statistical plan.
