#!/usr/bin/env python3
"""Verify start/stop codons for EVERY CDS in a PhageMine GenBank export.

Usage:
  python verify_codons.py genbank_submission/genome.gbk \
      --tbl genbank_submission/features.tbl \
      --val genbank_submission/table2asn_*/genome.val \
      --require-table 11

Exit code 0 only if every CDS passes every check. Exit code 1 otherwise.
Independent of PhageMine's own validator on purpose.
"""
import argparse
import glob
import json
import re
import sys

from Bio import SeqIO
from Bio.SeqFeature import AfterPosition, BeforePosition
from Bio.Seq import Seq

STOPS = {"TAA", "TAG", "TGA"}  # identical for NCBI tables 1 and 11
STARTS = {
    1: {"ATG", "CTG", "TTG"},
    11: {"TTG", "CTG", "ATT", "ATC", "ATA", "GTG", "ATG"},
}
# table2asn finding codes that mean a start/stop/partial problem
VAL_PATTERNS = [
    r"SEQ_FEAT\.StartCodon",
    r"SEQ_INST\.BadProteinStart",
    r"SEQ_FEAT\.NoStop",
    r"SEQ_FEAT\.InternalStop",
    r"SEQ_FEAT\.Partial\w*",
    r"SEQ_FEAT\.CDSwithNoStop\w*",
]


def partial_flags(loc, strand):
    """Return (partial5, partial3, left_partial, right_partial)."""
    left = isinstance(loc.start, BeforePosition)
    right = isinstance(loc.end, AfterPosition)
    if strand == -1:
        return right, left, left, right
    return left, right, left, right


def check_cds(rec, f, require_table, extended_starts=False):
    q = f.qualifiers
    tag = q.get("locus_tag", q.get("protein_id", ["?"]))[0]
    strand = f.location.strand or 1
    p5, p3, left, right = partial_flags(f.location, strand)
    problems = []

    declared = q.get("transl_table", [None])[0]
    table = int(declared) if declared else 1  # NCBI default is code 1
    if declared is None:
        problems.append("no transl_table declared (NCBI defaults to 1)")
    if require_table and table != require_table:
        problems.append(f"transl_table {table} != required {require_table}")

    cs = int(q.get("codon_start", ["1"])[0])
    if cs not in (1, 2, 3):
        problems.append(f"codon_start {cs} is not 1, 2 or 3")
    if not p5 and cs != 1:
        problems.append(f"complete 5-prime CDS has codon_start {cs}, expected 1")
    seq = f.extract(rec.seq)
    coding = seq[cs - 1:]
    first = str(coding[:3]).upper()
    last = str(coding[-3:]).upper()

    if len(coding) % 3 and not p3:
        problems.append(f"length {len(coding)} not a multiple of 3 and not 3' partial")

    # Start codon
    start_ok_set = STARTS.get(table, STARTS[1])
    if table == 11 and not extended_starts:
        start_ok_set = {"ATG", "GTG", "TTG"}
    if not p5 and first not in start_ok_set:
        problems.append(f"start {first} not valid for table {table} and CDS not 5' partial")

    # Stop codon
    if not p3 and last not in STOPS:
        problems.append(f"stop {last} is not TAA/TAG/TGA and CDS not 3' partial")

    # Partial markers only allowed at sequence ends
    if left and f.location.start != 0:
        problems.append(f"'<' at {int(f.location.start) + 1} but not at sequence start")
    if right and f.location.end != len(rec.seq):
        problems.append(f"'>' at {int(f.location.end)} but not at sequence end ({len(rec.seq)})")

    # Translation: no internal stops; qualifier matches
    usable = coding[: len(coding) - (len(coding) % 3)]
    prot = str(Seq(str(usable)).translate(table=table))
    body = prot[:-1] if (not p3 and prot.endswith("*")) else prot
    if "*" in body:
        problems.append(f"internal stop codon at aa {body.index('*') + 1}")
    if not p5 and prot and prot[0] != "M" and first in start_ok_set:
        # NCBI converts a valid alternative start to M; Biopython already does for table 11
        prot = "M" + prot[1:]
    given = q.get("translation", [None])[0]
    if given:
        expect = prot[:-1] if prot.endswith("*") else prot
        if given.upper() != expect.upper():
            problems.append("translation qualifier does not match the CDS sequence")

    else:
        problems.append("missing translation qualifier")

    return {
        "locus_tag": tag,
        "location": str(f.location),
        "strand": strand,
        "table": table,
        "start": first,
        "stop": last,
        "partial5": p5,
        "partial3": p3,
        "problems": problems,
    }


def check_tbl(path, require_table):
    """Every CDS block in features.tbl must carry transl_table."""
    problems, cds, tt = [], 0, 0
    in_cds = False
    for line in open(path):
        if line.strip() and not line.startswith("\t") and not line.startswith(">"):
            parts = line.rstrip("\n").split("\t")
            in_cds = len(parts) >= 3 and parts[2] == "CDS"
            cds += in_cds
        elif in_cds and re.match(r"^\t\t\ttransl_table\t", line):
            tt += 1
            if require_table and line.split("\t")[-1].strip() != str(require_table):
                problems.append(f"features.tbl has {line.strip()}")
    if cds != tt:
        problems.append(f"features.tbl: {cds} CDS rows but {tt} transl_table qualifiers")
    return problems


def check_val(paths):
    problems = []
    pat = re.compile("|".join(VAL_PATTERNS))
    for p in paths:
        for line in open(p):
            if pat.search(line):
                problems.append(f"{p.split('/')[-1]}: {line.strip()[:200]}")
    return problems


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gbk")
    ap.add_argument("--tbl")
    ap.add_argument("--val", nargs="*", default=[])
    ap.add_argument("--require-table", type=int, default=11)
    ap.add_argument("--json", help="write full results here")
    ap.add_argument("--extended-starts", action="store_true",
                    help="Accept all seven table-11 initiators; default accepts ATG/GTG/TTG")
    a = ap.parse_args()

    rec = SeqIO.read(a.gbk, "genbank")
    rows = [check_cds(rec, f, a.require_table, a.extended_starts) for f in rec.features if f.type == "CDS"]
    failed = [r for r in rows if r["problems"]]

    extra = []
    if a.tbl:
        extra += [f"[tbl] {p}" for p in check_tbl(a.tbl, a.require_table)]
    vals = [x for pat in a.val for x in glob.glob(pat)]
    extra += [f"[val] {p}" for p in check_val(vals)]

    print(f"Record {rec.id}: {len(rec.seq)} bp, {len(rows)} CDS checked")
    print(f"{'locus_tag':<12}{'location':<24}{'tbl':<4}{'start':<6}{'stop':<6}5'/3'  result")
    for r in rows:
        flag = f"{'P' if r['partial5'] else '-'}/{'P' if r['partial3'] else '-'}"
        res = "PASS" if not r["problems"] else "FAIL"
        if res == "FAIL" or a.json is None and len(rows) <= 15:
            print(f"{r['locus_tag']:<12}{r['location']:<24}{r['table']:<4}{r['start']:<6}{r['stop']:<6}{flag:<6}{res}")
        for p in r["problems"]:
            print(f"    - {p}")
    for e in extra:
        print(e)

    ok = not failed and not extra
    print(f"\n{len(rows) - len(failed)}/{len(rows)} CDS pass; {len(extra)} package-level problems")
    print("RESULT:", "PASS" if ok else "FAIL")
    if a.json:
        json.dump({"cds": rows, "package_problems": extra, "pass": ok}, open(a.json, "w"), indent=2)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
