# Public benchmark package

This compact package contains the frozen seven-phage benchmark summaries and
publication figures. It deliberately excludes raw tool outputs, annotation
logs, and third-party databases. The complete local evidence remains preserved
outside the public source package.

Tables use one-based inclusive coordinates and explicit strict/relaxed
structural definitions. Reference provenance records the accession-specific
GenBank inputs and checksums. Figures are derived from the TSV tables.

Verify the package with `shasum -a 256 -c SHA256SUMS` from the repository root.
