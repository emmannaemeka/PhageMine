# Scientific validation protocol

PhageMine cannot be declared validated by adding more product rules. A release
candidate intended for scientific comparison must pass this protocol.

## Locked benchmark

Use at least 100 expert-reviewed complete phage genomes spanning genome sizes,
hosts and major bacterial-virus families. Reserve incomplete contigs and
bacterial sequences as separate challenge sets. Freeze accession, sequence
checksum, annotation version, reviewer, review date and evidence for every
truth record.

Near-identical genomes and proteins present in PHROGs, VOGDB, Swiss-Prot or
INPHARED must be identified. Report both the full benchmark and a leakage-
controlled subset. Tune rules only on a development partition; report final
performance once on a held-out test partition.

## Comparators and metrics

Run Pharokka plus Phold, multiPhATE2, Prokka and PhageMine independently with
recorded versions and databases. Report exact CDS precision/recall/F1,
boundary-tolerant performance, product assertion precision, functional
coverage, abstention, hallmark sensitivity, unsafe-product rate, wall time and
peak memory. Stratify results by genome size and taxonomy.

Confidence calibration must be assessed separately for every evidence-strength
level. A label may be described as calibrated only after observed precision
and confidence intervals are published for held-out data.

## Reproducibility

Archive inputs, exact commands, environment lock, database manifests, raw
outputs, normalization vocabulary and analysis scripts. Tool agreement without
reviewed truth must be labelled `TOOL_AGREEMENT_ONLY`.
