# Database policy

The MVP ships no biological database. Its mock evidence adapter exists solely to demonstrate data flow and is visibly marked `MOCK` in all artifacts.

PhageMine may download databases with `phagemine databases install`, but the
downloaded data remain independent resources and are never embedded in the
software distribution. VOGDB and the PHROGs distribution are checksum-pinned.
Pfam and Swiss-Prot use provider current-release metadata and published
checksums, then record the resolved release and hashes in a local
`install_manifest.json`. A failed download, checksum, preparation step, or
readiness validation must not be registered as `READY`.

Production adapters should use versioned, locally documented snapshots of: curated protein records, broad similarity references, profile/domain libraries, and phage-focused protein-family collections. Every result must preserve database release, checksum, thresholds, coverage, identity, e-value, search tool/version, and provenance of its transferred annotation. Updates should be opt-in and create a new manifest because rankings can change.
