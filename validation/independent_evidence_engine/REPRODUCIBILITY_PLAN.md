# Reproducibility plan

The validation freeze manifest records repository commit and branch, tool and
provider versions, database versions and checksums, accession/version, FASTA
checksum, command line, configuration, thresholds, synonym-rule checksum,
random seeds, operating-system/Python environment, raw-output checksums,
reviewer-material checksums, analysis-script checksums, and generated-table
and figure checksums.

Raw predictions and evidence packets are immutable. Every transformation has a
deterministic script and a machine-readable schema. The tool-blinding key is
not included in public reviewer artifacts and is released only after explicit
unblinding authorization. Re-running from the manifest must fail loudly when a
required input, checksum, version, or adjudication field is missing.

The validation archive distinguishes development, validation, and external
test data. Validation genomes are never used as the sole tuning set. Reports
identify unavailable databases, failed providers, missing evidence, and
unresolved reference functions rather than substituting values.

`sample_size_scenarios.py` is the outcome-independent planning calculation. It
uses only declared variance, precision, effect-size, power, and confidence
assumptions and writes a JSON scenario table; it never reads a genome or tool
output.
