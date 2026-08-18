"""Deterministic, data-driven figures for annotation and discovery reports."""
from __future__ import annotations

import csv
import json
import os
from collections import Counter, defaultdict
from pathlib import Path


STATE_COLORS = {
    "KNOWN_FUNCTION": "#2166ac",
    "PROBABLE_FUNCTION": "#67a9cf",
    "FUNCTIONAL_CLASS_ONLY": "#d1e5f0",
    "CONSERVED_UNKNOWN": "#ef8a62",
    "CONFLICTING_EVIDENCE": "#b2182b",
    "UNRESOLVED": "#777777",
}


def _plotting(root: Path):
    os.environ.setdefault("MPLCONFIGDIR", str(root / ".matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def _save(fig, base: Path) -> list[str]:
    fig.savefig(base.with_suffix(".png"), dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    return [str(base.with_suffix(".png")), str(base.with_suffix(".svg"))]


def generate_annotation_figures(root, proteins, classifications, contexts=None, modules=None):
    root = Path(root); figures = root / "figures"; figures.mkdir(parents=True, exist_ok=True)
    data_dir = root / "figure_data"; data_dir.mkdir(exist_ok=True)
    manifest = {"created": [], "skipped": [], "failed": [], "source_data": []}
    if not proteins:
        manifest["skipped"].append({"figure": "all", "reason": "no predicted proteins"})
        (figures / "figure_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
        return manifest
    try:
        plt = _plotting(figures)
        states = Counter(row.get("functional_state") or "UNRESOLVED" for row in classifications)
        with (data_dir / "functional_classification_distribution.tsv").open("w") as h:
            h.write("functional_state\tcount\n"); h.writelines(f"{k}\t{v}\n" for k, v in sorted(states.items()))
        labels = [s for s in STATE_COLORS if states[s]]
        fig, ax = plt.subplots(figsize=(7.2, 4.5))
        ax.bar(labels, [states[s] for s in labels], color=[STATE_COLORS[s] for s in labels])
        ax.set(title="Functional classification of predicted proteins", ylabel="Proteins", xlabel="Functional state")
        ax.tick_params(axis="x", rotation=30); fig.tight_layout()
        manifest["created"] += _save(fig, figures / "functional_classification_distribution"); plt.close(fig)
        manifest["source_data"].append("functional_classification.tsv")

        sources = Counter(e.source for p in proteins for e in p.evidence if e.supports)
        with (data_dir / "evidence_support_matrix.tsv").open("w") as h:
            h.write("protein_id\tPfam\tVOGDB\tPHROGs\tSwiss-Prot\n")
            for p in proteins:
                present = {s: any(e.source == s and e.supports for e in p.evidence) for s in ("Pfam", "VOGDB", "PHROGs", "Swiss-Prot")}
                h.write(f"{p.protein_id}\t" + "\t".join("1" if present[s] else "0" for s in present) + "\n")
        if sources:
            fig, ax = plt.subplots(figsize=(6.5, 4.2)); names = sorted(sources)
            ax.bar(names, [sources[n] for n in names], color="#4c78a8")
            ax.set(title="Accepted evidence records by source", ylabel="Accepted records", xlabel="Evidence source")
            ax.tick_params(axis="x", rotation=25); fig.tight_layout()
            manifest["created"] += _save(fig, figures / "evidence_source_summary"); plt.close(fig)
            manifest["source_data"].append("evidence.json")
        else:
            manifest["skipped"].append({"figure": "evidence_source_summary", "reason": "no accepted evidence records"})

        by_id = {row["protein_id"]: row.get("functional_state") or "UNRESOLVED" for row in classifications}
        ordered = sorted(proteins, key=lambda p: (p.start, p.end, p.protein_id))
        fig, ax = plt.subplots(figsize=(12, max(2.5, min(5.5, len(ordered) / 25))))
        for p in ordered:
            direction = 1 if p.strand == "+" else -1
            x = p.start if direction == 1 else p.end
            dx = max(1, p.end - p.start) * direction
            ax.arrow(x, 0, dx, 0, width=.12, head_width=.34, head_length=max(8, abs(dx) * .18),
                     length_includes_head=True, color=STATE_COLORS.get(by_id.get(p.protein_id), "#999999"))
        ax.set(title="Genome functional map (PHANOTATE CDS boundaries preserved)", xlabel="Genome coordinate (nt)")
        ax.set_yticks([]); ax.set_ylim(-.7, .7); fig.tight_layout()
        manifest["created"] += _save(fig, figures / "genome_functional_map"); plt.close(fig)
        manifest["source_data"] += ["genes.gff3", "functional_classification.tsv"]
        manifest["source_data"] += ["figure_data/functional_classification_distribution.tsv", "figure_data/evidence_support_matrix.tsv"]
    except Exception as exc:
        manifest["failed"].append({"figure": "annotation figures", "reason": f"{type(exc).__name__}: {exc}"})
    (figures / "figure_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return manifest


def generate_discovery_figures(root, proteins, families, members, recurrence, ranking, neighbourhoods):
    root = Path(root); figures = root / "figures"; contexts_dir = figures / "genomic_context"
    figures.mkdir(parents=True, exist_ok=True); contexts_dir.mkdir(exist_ok=True)
    manifest = {"created": [], "skipped": [], "failed": [], "source_data": []}
    try:
        plt = _plotting(figures)
        states = Counter(p["functional_state"] for p in proteins)
        labels = [s for s in STATE_COLORS if states[s]]
        fig, ax = plt.subplots(figsize=(7.5, 4.8)); values = [states[s] for s in labels]
        bars = ax.bar(labels, values, color=[STATE_COLORS[s] for s in labels])
        total = sum(values)
        for bar, value in zip(bars, values): ax.text(bar.get_x()+bar.get_width()/2, value, f"{value}\n({value/total:.1%})", ha="center", va="bottom", fontsize=8)
        ax.set(title="Cohort functional-classification distribution", ylabel="Protein occurrences", xlabel="Functional state")
        ax.tick_params(axis="x", rotation=30); fig.tight_layout()
        manifest["created"] += _save(fig, figures / "functional_classification_distribution"); plt.close(fig)
        manifest["source_data"].append("cohort_unknown_proteome.tsv (targets); per-genome functional_classification.tsv (complete counts)")

        genomes = sorted({p["genome_id"] for p in proteins})
        recurrent = [r for r in recurrence if int(r["genome_count"]) > 1]
        heat_rows = recurrent[:80] if recurrent else recurrence[:80]
        if heat_rows and genomes:
            matrix = [[int(r["copies_by_genome"].get(g, 0)) for g in genomes] for r in heat_rows]
            height = max(4, min(18, len(heat_rows) * .24))
            fig, ax = plt.subplots(figsize=(max(7, len(genomes)*1.1), height)); image = ax.imshow(matrix, aspect="auto", cmap="Blues", interpolation="nearest")
            ax.set_xticks(range(len(genomes)), genomes, rotation=35, ha="right"); ax.set_yticks(range(len(heat_rows)), [r["pmf_id"] for r in heat_rows], fontsize=7)
            ax.set(title="PMF recurrence across genomes", xlabel="Genome", ylabel="PhageMine family (PMF)")
            fig.colorbar(image, ax=ax, label="Protein copies"); fig.tight_layout()
            manifest["created"] += _save(fig, figures / "pmf_recurrence_heatmap"); plt.close(fig)
            manifest["source_data"].append("family_recurrence.tsv")
        else:
            manifest["skipped"].append({"figure": "pmf_recurrence_heatmap", "reason": "no PMF-by-genome data"})

        distribution = Counter(int(r["genome_count"]) for r in recurrence)
        if distribution:
            xs = sorted(distribution); fig, ax = plt.subplots(figsize=(6.8, 4.5))
            ax.bar(xs, [distribution[x] for x in xs], color="#4c78a8")
            ax.set(title="PMF recurrence distribution", xlabel="Genomes represented", ylabel="PMFs", xticks=xs); fig.tight_layout()
            manifest["created"] += _save(fig, figures / "pmf_recurrence_distribution"); plt.close(fig)
            manifest["source_data"].append("family_recurrence.tsv")
        else:
            manifest["skipped"].append({"figure": "pmf_recurrence_distribution", "reason": "no PMFs"})

        top = ranking[:20]
        if top:
            fig, ax = plt.subplots(figsize=(8.5, max(4.5, len(top)*.33)))
            ys = list(range(len(top))); ax.barh(ys, [int(r["genome_count"]) for r in top], color="#e07a5f")
            ax.set_yticks(ys, [r["pmf_id"] for r in top]); ax.invert_yaxis()
            ax.set(title="Highest-priority recurrent unknown/unresolved PMFs", xlabel="Genomes represented", ylabel="PMF")
            fig.tight_layout(); manifest["created"] += _save(fig, figures / "discovery_ranking"); plt.close(fig)
            manifest["source_data"].append("discovery_ranking.tsv")
        else:
            manifest["skipped"].append({"figure": "discovery_ranking", "reason": "no recurrent unknown/unresolved PMFs met the existing priority rules"})

        member_by_family = defaultdict(list)
        for row in neighbourhoods: member_by_family[row["pmf_id"]].append(row)
        plotted = 0
        for candidate in ranking[:10]:
            rows = member_by_family.get(candidate["pmf_id"], [])
            if len({r["genome_id"] for r in rows}) < 2: continue
            fig, axes = plt.subplots(len(rows), 1, figsize=(11, max(2.5, len(rows)*1.25)), squeeze=False)
            for ax, row in zip(axes[:, 0], rows):
                neighborhood = row.get("neighbourhood", [])
                if not neighborhood: continue
                origin = min(n["start"] for n in neighborhood)
                for n in neighborhood:
                    direction = 1 if n["strand"] == "+" else -1; x = n["start"]-origin if direction == 1 else n["end"]-origin
                    dx = max(1, n["end"]-n["start"])*direction
                    color = "#f4a261" if n.get("is_target") else STATE_COLORS.get(n.get("functional_state"), "#999999")
                    ax.arrow(x, 0, dx, 0, width=.11, head_width=.3, head_length=max(5, abs(dx)*.15), length_includes_head=True, color=color)
                    ax.text((n["start"]+n["end"])/2-origin, .36, n.get("pmf_id") or n["protein_id"], ha="center", fontsize=6, rotation=25)
                ax.set_yticks([]); ax.set_ylabel(row["genome_id"], rotation=0, ha="right", va="center"); ax.set_ylim(-.55, .85)
            axes[-1,0].set_xlabel("Local genomic coordinate (nt; independently centered by neighborhood)")
            fig.suptitle(f"Gene neighborhoods for {candidate['pmf_id']} (context is supporting evidence)")
            fig.tight_layout(); manifest["created"] += _save(fig, contexts_dir / candidate["pmf_id"]); plt.close(fig); plotted += 1
        if plotted:
            manifest["source_data"].append("conserved_neighbourhoods.tsv")
        else:
            manifest["skipped"].append({"figure": "genomic_context/", "reason": "no prioritized PMF occurred in at least two genomes"})
        manifest["skipped"].append({"figure": "pmf_network", "reason": "no explicit, interpretable inter-family relationship model is produced; a network would imply unsupported edges"})
    except Exception as exc:
        manifest["failed"].append({"figure": "discovery figures", "reason": f"{type(exc).__name__}: {exc}"})
    (figures / "figure_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True))
    return manifest
