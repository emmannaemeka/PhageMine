# Database policy

The MVP ships no biological database. Its mock evidence adapter exists solely to demonstrate data flow and is visibly marked `MOCK` in all artifacts.

Production adapters should use versioned, locally documented snapshots of: curated protein records, broad similarity references, profile/domain libraries, and phage-focused protein-family collections. Every result must preserve database release, checksum, thresholds, coverage, identity, e-value, search tool/version, and provenance of its transferred annotation. Updates should be opt-in and create a new manifest because rankings can change.
