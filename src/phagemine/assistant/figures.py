"""Deterministic comparison figure: linear gene maps with ortholog links.

Every figure is written with its source table, so the picture can be checked
against data rather than trusted on sight.
"""
from __future__ import annotations

import csv
import os
from pathlib import Path

SHARED = "#2166ac"
UNIQUE = "#bdbdbd"
LINK = "#67a9cf"


def comparison_map(genomes: list[dict], links: list[dict], output_base: Path, title: str) -> list[Path]:
    """genomes: [{sample, length, genes:[{protein_id,start,end,strand,shared}]}]
    links: [{sample_a, protein_a, sample_b, protein_b}] between adjacent tracks."""
    output_base.parent.mkdir(parents=True, exist_ok=True)
    table = output_base.with_suffix(".tsv")
    with table.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["sample", "protein_id", "start", "end", "strand", "shared_with_another_genome"])
        for genome in genomes:
            for gene in genome["genes"]:
                writer.writerow([genome["sample"], gene["protein_id"], gene["start"], gene["end"], gene["strand"], gene["shared"]])
    os.environ.setdefault("MPLCONFIGDIR", str(output_base.parent / ".matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch, Polygon

    longest = max(g["length"] for g in genomes) or 1
    rows = len(genomes)
    fig, ax = plt.subplots(figsize=(10, 1.6 * rows + 1.2))
    track_y = {g["sample"]: (rows - i - 1) * 2.0 for i, g in enumerate(genomes)}
    positions = {}
    for genome in genomes:
        y = track_y[genome["sample"]]
        ax.plot([0, genome["length"]], [y, y], color="#555555", linewidth=1)
        ax.text(-longest * 0.01, y, genome["sample"], ha="right", va="center", fontsize=9)
        for gene in genome["genes"]:
            start, end = sorted((int(gene["start"]), int(gene["end"])))
            height = 0.35
            ax.add_patch(plt.Rectangle((start, y - height / 2 if gene["strand"] == "+" else y - height),
                                       end - start, height, facecolor=SHARED if gene["shared"] else UNIQUE,
                                       edgecolor="white", linewidth=1.2))
            positions[(genome["sample"], gene["protein_id"])] = ((start + end) / 2, y)
    for link in links:
        a = positions.get((link["sample_a"], link["protein_a"]))
        b = positions.get((link["sample_b"], link["protein_b"]))
        if a and b:
            ya, yb = sorted((a[1], b[1]), reverse=True)
            xa, xb = (a[0], b[0]) if a[1] >= b[1] else (b[0], a[0])
            ax.add_patch(Polygon([(xa, ya - 0.4), (xb, yb + 0.4)], closed=False, edgecolor=LINK, alpha=0.5, linewidth=0.8))
    ax.set_xlim(-longest * 0.02, longest * 1.02)
    ax.set_ylim(-1, (rows - 1) * 2.0 + 1)
    ax.set_yticks([])
    ax.set_xlabel("Genome position (bp)")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.set_title(title, fontsize=11, loc="left")
    ax.legend(handles=[Patch(color=SHARED, label="Homolog found in another genome"), Patch(color=UNIQUE, label="No homolog found"),
                       Patch(color=LINK, label="Homolog link")], loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=3, frameon=False, fontsize=8)
    fig.tight_layout()
    paths = [output_base.with_suffix(".png"), output_base.with_suffix(".svg")]
    fig.savefig(paths[0], dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(paths[1], bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return [*paths, table]
