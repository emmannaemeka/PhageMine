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
