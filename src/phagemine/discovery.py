"""Cohort-level biological discovery built from completed protein annotations."""
from __future__ import annotations

import csv
import hashlib
import html
import json
import os
import shutil
import subprocess
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from .family import build_database
from .family_external import validate_external
from .figures import generate_discovery_figures
from .inphared import compare_genomes


FAMILY_CONFIG = {"family_build": {"minimum_identity": 0.3, "minimum_coverage": 0.5,
                                   "coverage_mode": 0, "clustering_mode": 0}}
UNKNOWN_STATES = {"CONSERVED_UNKNOWN", "UNRESOLVED"}


def _sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""): digest.update(block)
    return digest.hexdigest()


def _signature(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), default=str).encode()).hexdigest()


def _read_tsv(path):
    with Path(path).open() as handle: return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv(path, rows, columns):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader(); writer.writerows(rows)


def _mmseqs_version(executable):
    exe = str(executable) if executable else shutil.which("mmseqs")
    if not exe or not (Path(exe).is_file() or shutil.which(exe)): return "unavailable"
    result = subprocess.run([exe, "version"], capture_output=True, text=True, check=False)
    return (result.stdout or result.stderr).strip().splitlines()[0] or "unknown"


def _load_samples(sample_dirs):
    proteins = []; samples = []; fingerprints = {}
    required = ("evidence.json", "functional_classification.json", "genomic_context.json", "run_manifest.json", "proteins.faa")
    for sample in sorted(map(Path, sample_dirs), key=lambda p: p.name):
        missing = [name for name in required if not (sample / name).is_file()]
        if missing: raise ValueError(f"Discovery input {sample} is incomplete: {', '.join(missing)}")
        fingerprints[str(sample)] = {name: _sha(sample / name) for name in required}
        evidence = json.loads((sample / "evidence.json").read_text())
        classifications = {r["protein_id"]: r for r in json.loads((sample / "functional_classification.json").read_text())}
        contexts = {r["protein_id"]: r for r in json.loads((sample / "genomic_context.json").read_text())}
        manifest = json.loads((sample / "run_manifest.json").read_text())
        genome_ids = sorted({p.get("genome_id") for p in evidence if p.get("genome_id")})
        if len(genome_ids) != 1: raise ValueError(f"Discovery input {sample} does not contain exactly one stable genome ID")
        genome_id = genome_ids[0]
        samples.append({"sample_id": sample.name, "genome_id": genome_id, "path": str(sample),
                        "input_sha256": manifest.get("input_sha256")})
        for p in evidence:
            local_id = p["protein_id"]; classification = classifications.get(local_id, {}); context = contexts.get(local_id, {})
            accepted = [e for e in p.get("evidence", []) if e.get("supports")]
            sources = sorted({e.get("source") for e in accepted if e.get("source")})
            summaries = sorted({f"{e.get('source')}:{e.get('identifier') or e.get('family_name') or e.get('description')}" for e in accepted})
            proteins.append({
                "cohort_protein_id": f"{genome_id}|{local_id}", "genome_id": genome_id,
                "sample_id": sample.name, "protein_id": local_id, "start": int(p["start"]), "end": int(p["end"]),
                "strand": p["strand"], "sequence": p["sequence"], "sequence_sha256": hashlib.sha256(p["sequence"].encode()).hexdigest(),
                "annotation": p.get("annotation"), "functional_state": classification.get("functional_state", "UNRESOLVED"),
                "proposed_function": classification.get("proposed_function"), "functional_category": classification.get("functional_category"),
                "confidence": classification.get("confidence"), "conservation_status": classification.get("conservation_status"),
                "evidence_sources": sources, "evidence_summary": summaries, "context": context,
            })
    identifiers = [p["cohort_protein_id"] for p in proteins]
    if len(identifiers) != len(set(identifiers)): raise ValueError("Cohort protein IDs are not unique after genome qualification")
    return proteins, samples, fingerprints


def _family_tables(families, raw_members, proteins):
    by_id = {p["cohort_protein_id"]: p for p in proteins}; grouped = defaultdict(list)
    for member in raw_members: grouped[member["family_id"]].append(member)
    family_rows = []; member_rows = []; recurrence = []
    for family in sorted(families, key=lambda x: x["family_id"]):
        fid = family["family_id"]; ms = grouped[fid]; ps = [by_id[m["member_id"]] for m in ms]
        states = Counter(p["functional_state"] for p in ps); genomes = sorted({p["genome_id"] for p in ps})
        copies = Counter(p["genome_id"] for p in ps); sources = sorted({s for p in ps for s in p["evidence_sources"]})
        evidence = sorted({s for p in ps for s in p["evidence_summary"]})
        characterized = sum(p["functional_state"] not in UNKNOWN_STATES for p in ps)
        unknown = len(ps) - characterized
        dominant = sorted(states, key=lambda state: (-states[state], state))[0]
        family_rows.append({"pmf_id": fid, "representative_protein": family["representative_protein_id"],
            "member_count": len(ps), "genome_count": len(genomes), "genomes_represented": ",".join(genomes),
            "member_protein_ids": ",".join(sorted(p["cohort_protein_id"] for p in ps)),
            "functional_state_composition": json.dumps(dict(sorted(states.items())), sort_keys=True),
            "dominant_functional_state": dominant, "characterized_member_count": characterized,
            "unknown_unresolved_member_count": unknown, "evidence_sources": ",".join(sources),
            "available_evidence_summary": ";".join(evidence)})
        recurrence.append({"pmf_id": fid, "member_count": len(ps), "genome_count": len(genomes),
            "genomes_represented": ",".join(genomes), "within_genome_duplicate_count": sum(max(0, n-1) for n in copies.values()),
            "copies_by_genome": dict(sorted(copies.items())), "recurrence_state": "RECURRENT_ACROSS_GENOMES" if len(genomes)>1 else ("REPEATED_WITHIN_ONE_GENOME" if len(ps)>1 else "SINGLE_GENOME_SINGLE_COPY")})
        for p in sorted(ps, key=lambda x: (x["genome_id"], x["start"], x["cohort_protein_id"])):
            member_rows.append({"pmf_id": fid, "cohort_protein_id": p["cohort_protein_id"], "genome_id": p["genome_id"],
                "protein_id": p["protein_id"], "start": p["start"], "end": p["end"], "strand": p["strand"],
                "sequence_sha256": p["sequence_sha256"], "functional_state": p["functional_state"],
                "proposed_function": p["proposed_function"], "evidence_sources": ",".join(p["evidence_sources"])})
    return family_rows, member_rows, recurrence


def _neighbourhoods(proteins, member_rows):
    by_global = {p["cohort_protein_id"]: p for p in proteins}; pmf_by_protein = {m["cohort_protein_id"]: m["pmf_id"] for m in member_rows}
    ordered = defaultdict(list)
    for p in proteins: ordered[p["genome_id"]].append(p)
    for genome in ordered: ordered[genome].sort(key=lambda x: (x["start"], x["end"], x["protein_id"]))
    index = {p["cohort_protein_id"]: i for genome in ordered for i,p in enumerate(ordered[genome])}
    occurrences = []; signatures = defaultdict(lambda: defaultdict(list))
    copies_seen = Counter()
    for member in member_rows:
        target = by_global[member["cohort_protein_id"]]; genome = target["genome_id"]; i = index[target["cohort_protein_id"]]
        neighborhood = []
        for j in range(max(0, i-2), min(len(ordered[genome]), i+3)):
            p = ordered[genome][j]
            neighborhood.append({"relative_gene_offset": j-i, "cohort_protein_id": p["cohort_protein_id"], "protein_id": p["protein_id"],
                "pmf_id": pmf_by_protein.get(p["cohort_protein_id"]), "start": p["start"], "end": p["end"], "strand": p["strand"],
                "functional_state": p["functional_state"], "proposed_function": p["proposed_function"], "functional_category": p["functional_category"],
                "is_target": j == i})
        signature={(n["relative_gene_offset"],n["pmf_id"]) for n in neighborhood if not n["is_target"] and n["pmf_id"]}
        signatures[member["pmf_id"]][genome].append(signature); copies_seen[(member["pmf_id"],genome)] += 1
        occurrences.append({"pmf_id": member["pmf_id"], "cohort_protein_id": target["cohort_protein_id"], "genome_id": genome,
            "protein_id": target["protein_id"], "copy_index_within_genome": copies_seen[(member["pmf_id"],genome)], "neighbourhood": neighborhood})
    status_by_family = {}
    for fid, genome_sets in signatures.items():
        collapsed = {genome: set().union(*sets) for genome,sets in genome_sets.items()}
        union = set().union(*collapsed.values()) if collapsed else set()
        common = set.intersection(*collapsed.values()) if collapsed else set()
        if len(collapsed)<2: status="INSUFFICIENT_CONTEXT"
        elif common and all(value == next(iter(collapsed.values())) for value in collapsed.values()): status="STRONGLY_CONSERVED_CONTEXT"
        elif common: status="CONSERVED_CONTEXT"
        else: status="NO_CONSERVED_CONTEXT"
        status_by_family[fid]={"context_conservation_state":status, "context_conservation_support":len(common)/len(union) if union else 0.0,
            "shared_offset_pmf_anchors": sorted(f"{offset}:{pmf}" for offset,pmf in common), "genomes_compared":len(collapsed)}
    for row in occurrences: row.update(status_by_family[row["pmf_id"]])
    return occurrences, status_by_family


def _ranking(families, context_status, validation):
    external = {r["family_id"]: r for r in validation}; candidates=[]
    for family in families:
        composition=json.loads(family["functional_state_composition"])
        if family["genome_count"] <= 1 or sum(composition.get(s,0) for s in UNKNOWN_STATES) != family["member_count"]: continue
        ext=external.get(family["pmf_id"], {}); ctx=context_status.get(family["pmf_id"], {})
        state = "CONSERVED_UNKNOWN" if composition.get("CONSERVED_UNKNOWN") else "UNRESOLVED"
        pmfdb = ext.get("external_validation_status", "PMFDB_UNAVAILABLE")
        rationale = f"recurrent unknown/unresolved family across {family['genome_count']} genomes"
        rationale += f"; genomic context {ctx.get('context_conservation_state','INSUFFICIENT_CONTEXT').lower()}"
        rationale += f"; PMFDB state {pmfdb}; database absence is not interpreted as novelty"
        candidates.append({"pmf_id":family["pmf_id"], "functional_state":state, "member_count":family["member_count"],
            "genome_count":family["genome_count"], "genomes_represented":family["genomes_represented"],
            "context_conservation_state":ctx.get("context_conservation_state","INSUFFICIENT_CONTEXT"),
            "context_conservation_support":ctx.get("context_conservation_support",0.0), "pmfdb_state":pmfdb,
            "characterized_homolog_count":int(ext.get("characterized_homolog_count") or 0),
            "evidence_sources":family["evidence_sources"], "discovery_priority":"PRIORITIZED_RECURRENT_UNKNOWN",
            "reason_for_priority":rationale})
    candidates.sort(key=lambda r:(-r["genome_count"],-r["member_count"],-len([x for x in r["evidence_sources"].split(',') if x]),r["pmf_id"]))
    for rank,row in enumerate(candidates,1): row["rank"]=rank
    return candidates


def _reports(root, samples, proteins, families, ranking, figures, evidence_reused, pmfdb_provenance, inphared_comparison):
    counts=Counter(p["functional_state"] for p in proteins)
    lines=["# PhageMine discovery report", "", "> PMF recurrence, genomic context, and database matches are computational evidence. No external match is not evidence of novelty.", "",
        f"- Genomes: **{len(samples)}**", f"- Protein occurrences: **{len(proteins)}**", f"- PMFs: **{len(families)}**",
        f"- Recurrent PMFs: **{sum(int(f['genome_count'])>1 for f in families)}**", f"- Evidence reused from Annotation Mode: **{evidence_reused}**", "", "## Functional states", ""]
    lines += [f"- {state}: {counts[state]}" for state in sorted(counts)]
    lines += ["", "## Highest-priority recurrent unknown/unresolved PMFs", "", "| Rank | PMF | State | Members | Genomes | Context | PMFDB | Reason |", "|---:|---|---|---:|---:|---|---|---|"]
    for row in ranking[:50]:
        lines.append(f"| {row['rank']} | {row['pmf_id']} | {row['functional_state']} | {row['member_count']} | {row['genome_count']} | {row['context_conservation_state']} | {row['pmfdb_state']} | {row['reason_for_priority']} |")
    lines += ["", "## Figures", ""]
    for path in figures.get("created",[]): lines.append(f"- [{Path(path).name}]({Path(path).relative_to(root)})")
    for item in figures.get("skipped",[]): lines.append(f"- Skipped `{item['figure']}`: {item['reason']}")
    if pmfdb_provenance: lines += ["", "## PMFDB provenance", "", f"- Version: `{pmfdb_provenance.get('pmfdb_version')}`", f"- Reference: `{pmfdb_provenance.get('reference')}`"]
    nearest=(inphared_comparison or {}).get("unique_reference_summaries") or (inphared_comparison or {}).get("matches") or []
    if nearest:
        lines += ["", "## Whole-genome numerical taxonomy: INPHARED nearest references", "", "> Mash selects candidates only. Similarity is calculated using a VIRIDIC-compatible bidirectional BLASTN method; threshold agreement is not a formal ICTV assignment.", "", "| Query | Rank | Reference | Intergenomic similarity (%) | Query aligned (%) | Reference aligned (%) | Host | Taxon | Interpretation |", "|---|---:|---|---:|---:|---:|---|---|---|"]
        for row in nearest:
            taxon=row.get("phage_genus") or row.get("phage_family") or ""
            reference=row.get("reference_description") or row.get("reference_accession") or ""
            distance=row.get("best_mash_distance", row.get("mash_distance"))
            rank=row.get("rank", "unique")
            similarity=row.get("intergenomic_similarity_percent")
            lines.append(f"| {row['sample_id']} | {rank} | {reference} | {similarity if similarity is not None else 'Not calculated'} | {row.get('query_aligned_percent') if row.get('query_aligned_percent') is not None else 'Not calculated'} | {row.get('reference_aligned_percent') if row.get('reference_aligned_percent') is not None else 'Not calculated'} | {row.get('host_genus') or ''} | {taxon} | {row.get('taxonomic_interpretation') or 'Mash screening only'} |")
    markdown="\n".join(lines)+"\n"; (root/"discovery_report.md").write_text(markdown)
    figure_html="".join(f'<figure><a href="{html.escape(str(Path(path).relative_to(root)))}"><img src="{html.escape(str(Path(path).relative_to(root)))}" style="max-width:100%"></a><figcaption>{html.escape(Path(path).stem)}</figcaption></figure>' for path in figures.get("created",[]) if path.endswith('.png'))
    table_rows="".join(f"<tr><td>{r['rank']}</td><td>{r['pmf_id']}</td><td>{r['functional_state']}</td><td>{r['member_count']}</td><td>{r['genome_count']}</td><td>{r['context_conservation_state']}</td><td>{r['pmfdb_state']}</td></tr>" for r in ranking[:50])
    (root/"discovery_report.html").write_text(f"<!doctype html><html><head><meta charset='utf-8'><title>PhageMine discovery report</title><style>body{{font-family:Arial,sans-serif;max-width:1200px;margin:auto;padding:2rem}}table{{border-collapse:collapse}}th,td{{border:1px solid #bbb;padding:.35rem}}figure{{margin:2rem 0}}</style></head><body><h1>PhageMine discovery report</h1><p><strong>Scientific caution:</strong> context and homology are supporting evidence; NO_EXTERNAL_MATCH does not mean novel.</p><p>Genomes: {len(samples)}; proteins: {len(proteins)}; PMFs: {len(families)}.</p><table><thead><tr><th>Rank</th><th>PMF</th><th>State</th><th>Members</th><th>Genomes</th><th>Context</th><th>PMFDB</th></tr></thead><tbody>{table_rows}</tbody></table><h2>Figures</h2>{figure_html}</body></html>")


def build_discovery_outputs(sample_dirs, output, *, mmseqs="mmseqs", pmfdb=None, inphared=None, mash="mash", progress=None,
                            resume_existing=False, evidence_reused=False, source_mode="discover"):
    root=Path(output); root.mkdir(parents=True, exist_ok=True); started=time.monotonic(); timings={}
    checkpoints_path=root/"discovery_checkpoints.json"
    checkpoints=json.loads(checkpoints_path.read_text()) if resume_existing and checkpoints_path.is_file() else {}
    def emit(state,name,started_at=None,detail=""):
        elapsed=time.monotonic()-started_at if started_at is not None else None; suffix=f" ({elapsed:.1f}s)" if elapsed is not None else ""; suffix += f"; {detail}" if detail else ""
        if progress: progress._write(f"[Discovery] {state}: {name}{suffix}")
        if elapsed is not None: timings[name]=elapsed
    def mark(name, signature, paths, reused=False):
        checkpoints[name]={"state":"COMPLETE", "signature":signature, "reused":reused,
            "time":datetime.now(timezone.utc).isoformat(), "output_checksums":{str(Path(p).relative_to(root)):_sha(p) for p in paths if Path(p).is_file()}}
        checkpoints_path.write_text(json.dumps(checkpoints,indent=2,sort_keys=True))
    proteins,samples,fingerprints=_load_samples(sample_dirs); input_signature=_signature(fingerprints)
    pooled_fasta=root/"pooled_proteins.faa"; pooled_fasta.write_text("".join(f">{p['cohort_protein_id']}\n{p['sequence']}\n" for p in proteins))
    metadata={p["cohort_protein_id"]:{"genome_id":p["genome_id"],"protein_id":p["protein_id"],"annotation":p["proposed_function"] or p["annotation"]} for p in proteins}

    stage="PMF clustering"; t=time.monotonic(); emit("RUNNING",stage)
    family_db=root/"pmf_database"; cluster_sig=_signature({"inputs":fingerprints,"config":FAMILY_CONFIG,"backend":"MMSEQS2","mmseqs":_mmseqs_version(mmseqs)})
    reuse=resume_existing and checkpoints.get(stage,{}).get("signature")==cluster_sig and (family_db/"families.json").is_file() and (family_db/"family_members.tsv").is_file()
    if reuse: families=json.loads((family_db/"families.json").read_text())
    else: families=build_database(pooled_fasta,family_db,metadata=metadata,database_version="1.0",backend="MMSEQS2",config=FAMILY_CONFIG,mmseqs=mmseqs)
    raw_members=_read_tsv(family_db/"family_members.tsv"); family_rows,member_rows,recurrence=_family_tables(families,raw_members,proteins)
    mark(stage,cluster_sig,[family_db/"families.json",family_db/"family_members.tsv"],reuse); emit("DONE",stage,t,"reused" if reuse else f"{len(family_rows)} PMFs")

    stage="Recurrence analysis"; t=time.monotonic(); emit("RUNNING",stage)
    family_cols=["pmf_id","representative_protein","member_count","genome_count","genomes_represented","member_protein_ids","functional_state_composition","dominant_functional_state","characterized_member_count","unknown_unresolved_member_count","evidence_sources","available_evidence_summary"]
    member_cols=["pmf_id","cohort_protein_id","genome_id","protein_id","start","end","strand","sequence_sha256","functional_state","proposed_function","evidence_sources"]
    recurrence_cols=["pmf_id","member_count","genome_count","genomes_represented","within_genome_duplicate_count","copies_by_genome","recurrence_state"]
    _write_tsv(root/"pmf_families.tsv",family_rows,family_cols); _write_tsv(root/"pmf_members.tsv",member_rows,member_cols)
    recurrence_tsv=[{**r,"copies_by_genome":json.dumps(r["copies_by_genome"],sort_keys=True)} for r in recurrence]
    _write_tsv(root/"family_recurrence.tsv",recurrence_tsv,recurrence_cols)
    (root/"pmf_families.json").write_text(json.dumps(family_rows,indent=2,sort_keys=True)); (root/"pmf_members.json").write_text(json.dumps(member_rows,indent=2,sort_keys=True)); (root/"family_recurrence.json").write_text(json.dumps(recurrence,indent=2,sort_keys=True))
    unknown=[p for p in proteins if p["functional_state"] in UNKNOWN_STATES]
    unknown_cols=["cohort_protein_id","genome_id","protein_id","start","end","strand","sequence_sha256","functional_state","proposed_function","functional_category","confidence","conservation_status","evidence_sources","evidence_summary"]
    _write_tsv(root/"cohort_unknown_proteome.tsv",[{**p,"evidence_sources":",".join(p["evidence_sources"]),"evidence_summary":";".join(p["evidence_summary"])} for p in unknown],unknown_cols)
    mark(stage,input_signature,[root/"pmf_families.tsv",root/"pmf_members.tsv",root/"family_recurrence.tsv"]); emit("DONE",stage,t,f"{sum(r['genome_count']>1 for r in recurrence)} recurrent PMFs")

    stage="Genomic-context analysis"; t=time.monotonic(); emit("RUNNING",stage)
    neighbourhoods,context_status=_neighbourhoods(proteins,member_rows)
    neighbourhood_cols=["pmf_id","cohort_protein_id","genome_id","protein_id","copy_index_within_genome","genomes_compared","context_conservation_state","context_conservation_support","shared_offset_pmf_anchors","neighbourhood"]
    neighbourhood_tsv=[{**r,"shared_offset_pmf_anchors":",".join(r["shared_offset_pmf_anchors"]),"neighbourhood":json.dumps(r["neighbourhood"],sort_keys=True)} for r in neighbourhoods]
    _write_tsv(root/"conserved_neighbourhoods.tsv",neighbourhood_tsv,neighbourhood_cols); (root/"conserved_neighbourhoods.json").write_text(json.dumps(neighbourhoods,indent=2,sort_keys=True))
    mark(stage,input_signature,[root/"conserved_neighbourhoods.tsv"]); emit("DONE",stage,t)

    stage="PMFDB validation"; t=time.monotonic(); emit("RUNNING",stage)
    pmfdb_provenance={}; validation=[]
    pmf_sig=_signature({"families":[(f["pmf_id"],f["representative_protein"]) for f in family_rows],"pmfdb":str(pmfdb),"mmseqs":_mmseqs_version(mmseqs)})
    reuse_pmf=resume_existing and checkpoints.get(stage,{}).get("signature")==pmf_sig and (root/"pmfdb_validation.tsv").is_file()
    if reuse_pmf:
        validation=_read_tsv(root/"pmfdb_validation.tsv")
        if (root/"pmfdb_validation.json").is_file(): pmfdb_provenance=json.loads((root/"pmfdb_validation.json").read_text()).get("provenance",{})
    elif pmfdb:
        validation=validate_external(family_db,output=root/"pmfdb",mmseqs=mmseqs,reference_db=pmfdb)
        payload=json.loads((root/"pmfdb"/"pmf_external_validation.json").read_text()); pmfdb_provenance=payload["provenance"]
        shutil.copyfile(root/"pmfdb"/"pmf_external_validation.tsv",root/"pmfdb_validation.tsv")
        (root/"pmfdb_validation.json").write_text(json.dumps({"provenance":pmfdb_provenance,"families":validation},indent=2,sort_keys=True))
    else:
        validation=[{"family_id":f["pmf_id"],"query_representative":f["representative_protein"],"external_match_count":0,"characterized_homolog_count":0,"uncharacterized_homolog_count":0,"external_validation_status":"PMFDB_UNAVAILABLE"} for f in family_rows]
        _write_tsv(root/"pmfdb_validation.tsv",validation,list(validation[0]) if validation else ["family_id"])
        (root/"pmfdb_validation.json").write_text(json.dumps({"provenance":{"status":"PMFDB_UNAVAILABLE"},"families":validation},indent=2,sort_keys=True))
    mark(stage,pmf_sig,[root/"pmfdb_validation.tsv",root/"pmfdb_validation.json"],reuse_pmf); emit("DONE" if pmfdb or reuse_pmf else "SKIPPED",stage,t,"reused" if reuse_pmf else (pmfdb_provenance.get("pmfdb_version") or "PMFDB not configured"))

    stage="INPHARED genome comparison"; t=time.monotonic(); emit("RUNNING",stage)
    inphared_comparison={"status":"INPHARED_UNAVAILABLE","matches":[]}
    if inphared:
        provenance=inphared.get("provenance") or {}
        fastas={Path(sample).name:Path(sample)/"analysis_genome.fasta" for sample in sample_dirs if (Path(sample)/"analysis_genome.fasta").is_file()}
        inphared_comparison=compare_genomes(
            fastas,
            mash_index=provenance.get("mash_index_path"),
            metadata=provenance.get("metadata_path"),
            output=root,
            mash=mash,
            reference_fasta=inphared.get("path"),
        )
        inphared_comparison["resource_version"]=inphared.get("version")
        inphared_comparison["resource_manifest"]=provenance.get("reference_manifest_path")
        (root/"inphared_nearest_phages.json").write_text(json.dumps(inphared_comparison,indent=2,sort_keys=True)+"\n")
        emit("DONE",stage,t,f"{len(inphared_comparison.get('matches') or [])} nearest-reference rows")
    else:
        _write_tsv(root/"inphared_nearest_phages.tsv",[],["sample_id","rank","reference_accession","mash_distance","p_value","matching_hashes","reference_description","host_genus","phage_genus","phage_subfamily","phage_family","intergenomic_similarity_percent","query_aligned_percent","reference_aligned_percent","genome_length_ratio","similarity_method","species_threshold_percent","genus_threshold_percent","threshold_source","taxonomic_interpretation","interpretation"])
        _write_tsv(root/"inphared_summary.tsv",[],["sample_id","relationship","reference_description","reference_accessions","best_mash_distance","matching_hashes","intergenomic_similarity_percent","query_aligned_percent","reference_aligned_percent","genome_length_ratio","similarity_method","species_threshold_percent","genus_threshold_percent","threshold_source","taxonomic_interpretation","interpretation"])
        (root/"inphared_nearest_phages.json").write_text(json.dumps(inphared_comparison,indent=2,sort_keys=True)+"\n")
        emit("SKIPPED",stage,t,"INPHARED genomes not configured")
    mark(stage,_signature({"inputs":fingerprints,"inphared":str((inphared or {}).get('path'))}),[root/"inphared_nearest_phages.tsv",root/"inphared_summary.tsv",root/"inphared_nearest_phages.json"])

    stage="Discovery ranking"; t=time.monotonic(); emit("RUNNING",stage)
    ranking=_ranking(family_rows,context_status,validation)
    rank_cols=["rank","pmf_id","functional_state","member_count","genome_count","genomes_represented","context_conservation_state","context_conservation_support","pmfdb_state","characterized_homolog_count","evidence_sources","discovery_priority","reason_for_priority"]
    _write_tsv(root/"discovery_ranking.tsv",ranking,rank_cols); (root/"discovery_ranking.json").write_text(json.dumps(ranking,indent=2,sort_keys=True))
    mark(stage,input_signature,[root/"discovery_ranking.tsv"]); emit("DONE",stage,t,f"{len(ranking)} prioritized recurrent unknown/unresolved PMFs")

    stage="Figure generation"; t=time.monotonic(); emit("RUNNING",stage)
    figures=generate_discovery_figures(root,proteins,family_rows,member_rows,recurrence,ranking,neighbourhoods)
    emit("DONE" if not figures["failed"] else "FAILED",stage,t,f"{len(figures['created'])} files created")
    stage="Final output generation"; t=time.monotonic(); emit("RUNNING",stage)
    _reports(root,samples,proteins,family_rows,ranking,figures,evidence_reused,pmfdb_provenance,inphared_comparison)
    manifest={"pipeline":"PhageMine","mode":source_mode,"created_at":datetime.now(timezone.utc).isoformat(),"input_fingerprints":fingerprints,
        "samples":samples,"genome_count":len(samples),"total_proteins":len(proteins),"exact_unique_proteins":len({p['sequence_sha256'] for p in proteins}),
        "pmf_count":len(family_rows),"recurrent_pmf_count":sum(r["genome_count"]>1 for r in recurrence),"evidence_reused":evidence_reused,
        "family_clustering":{"backend":"MMSEQS2","thresholds":FAMILY_CONFIG["family_build"],"mmseqs_version":_mmseqs_version(mmseqs)},
        "pmfdb":pmfdb_provenance or {"status":"PMFDB_UNAVAILABLE"},"inphared":inphared_comparison,"figures":figures,"stage_timings_seconds":timings,
        "scientific_interpretation":{"no_external_match":"does not mean novel","genomic_context":"supporting evidence, not proof of function"},
        "total_runtime_seconds":time.monotonic()-started}
    (root/"pooled_manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)); mark(stage,input_signature,[root/"pooled_manifest.json",root/"discovery_report.html",root/"discovery_report.md"]); emit("DONE",stage,t)
    return manifest
