# Installation

Install the release candidate in a clean environment:

```bash
python -m venv .venv
. .venv/bin/activate
python -m pip install phagemine-1.0.0-py3-none-any.whl
```

Production gene prediction requires the external PHANOTATE script. Pass its
path explicitly when it is not available as `phanotate` on `PATH`, for example
`--phanotate /path/to/phanotate.py`.

Minimal Core run:

```bash
phagemine run genome.fasta --output results/run \
  --phanotate /path/to/phanotate.py
```

Optional evidence executables and databases are installed separately and are
not required for Core. Missing optional resources are reported and the Core
workflow remains usable.
