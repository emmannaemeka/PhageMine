# Automated evidence adjudication

This directory is a separate `AUTOMATED_EVIDENCE_ADJUDICATION` analysis track
for the frozen v1.2 benchmark. It does not alter the original human-review
protocol, blank workbook, key, predictions, reference records, matching rules,
or utility definition. Decisions are deterministic and use only the neutral
frozen evidence packets. The output must not be described as expert human
adjudication.

Run `python scripts/run_automated.py blind` to create and checksum the blind
outputs. Only after that freeze is verified may `python scripts/run_automated.py
unblind` read the confidential mapping and generate descriptive tool-level
summaries. Unresolvable cases remain missing; they are never silently scored as
incorrect or zero.
