# Independent panel-selection protocol

Panel selection occurs after this protocol is reviewed and frozen. Candidate
metadata may be screened for eligibility, but no candidate is run through
PhageMine during selection.

## Required exclusions

Exclude the five development genomes, all seven v1.2 reference genomes,
genomes used to design K2/K3/K4, genomes used for threshold tuning or
debugging, and close duplicates of any excluded genome. Exclusion decisions
are recorded with accession/version, sequence checksum, and provenance.

## Inclusion targets

Prefer complete, well-curated bacteriophage genomes with an independently
defensible reference record. Seek balanced representation of tailed dsDNA,
filamentous, small ssDNA, temperate/integration-capable, lytic, compact,
large, overlap-containing, short-CDS-containing, and nucleotide-divergent
phages. No single family, host, curator, or annotation source should dominate.

The panel is stratified prospectively by architecture, genome size, reference
quality, and database-overlap status. A candidate is not retained merely
because it is difficult or likely to favour one tool.

## Practical target and minimum representation

The planning target is a minimum analyzable panel of 32 independent genomes,
a preferred target of 64, and a reserve pool of 16 eligible genomes. These are
planning targets, not observed results: `sample_size_scenarios.py` evaluates
precision and paired-effect scenarios without validation data. The minimum
supports genome-level exploratory estimation and a paired comparison under
moderate planning variance; the preferred target also supports useful
precision under a conservative variance scenario. If exclusions or clustering
reduce the analyzable count, the reserve pool is used before any outcome is
inspected.

For architecture, seek 8–16 genomes in each common supported architecture and
at least four in each represented non-rare architecture stratum. Equal counts
are not forced for rare architectures; failure to reach four is reported as a
design limitation. For genome size, seek at least 20% compact (<50 kb), 20%
large (>150 kb), and the remainder intermediate (50–150 kb), subject to
curation availability. Seek at least 25% lytic and 25% temperate/
integration-capable genomes, plus at least four genomes with validated short or
overlapping CDSs and at least four with no significant nucleotide reference.
These are targets and minimums, not quotas that justify including weakly
curated genomes.

## Partitioning

If any development/tuning data are used to construct analysis code, they are
kept in a development partition. The validation test partition is untouched
until the protocol, software commit, databases, panel, and analysis scripts
are frozen. No validation outcome may influence inclusion, thresholds, or
synonym rules.

## Metadata to capture

Accession/version, source, completeness, architecture evidence, host/context,
sequence checksum, duplicate/cluster identifiers, reference quality,
database-overlap flags, exclusion rationale, and panel stratum are recorded
in a signed panel manifest. The manifest contains no PhageMine output before
the validation freeze.

## Relatedness control

Cluster candidate genomes before selection using nucleotide ANI **and** aligned
coverage jointly: genomes are in the same primary cluster when ANI is ≥95% and
aligned coverage is ≥80%. One representative is selected deterministically by
complete sequence, highest reference-standard tier, then stable accession
order; tool outputs and expected difficulty are never selection criteria. A
sensitivity clustering rule uses ANI ≥90% and coverage ≥70%. Cluster IDs,
pairwise values, methods, and excluded near-duplicates are retained in the
manifest. Protein-family clustering is a secondary leakage descriptor and does
not replace genome clustering.

## Comparator role

Pharokka may be run as an established comparator for paired structural and
functional outcomes defined for both tools. It is not penalized for
PhageMine-only evidence objects. PhageMine confidence, modules, architecture,
comparative significance, provenance, and review triage are evaluated
independently; any paired software comparison remains descriptive unless a
separate estimand is approved before execution.
