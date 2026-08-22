#!/bin/sh

cat <<'EOF'

PhageMine was installed successfully.

The biological evidence databases are distributed separately. Install and
validate them before running full evidence mode:

  phagemine databases install --all
  phagemine doctor

Individual installers are available for pfam, vogdb, swissprot, phrogs,
pmfdb and inphared.

EOF
