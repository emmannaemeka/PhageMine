"""PhageMine GUI — local researcher interface over the native CLI artifacts."""
from __future__ import annotations

import os
import threading
from pathlib import Path

import streamlit as st

from phagemine.gui import GUI_VERSION
from phagemine.gui.services.execution import build_cli_args, display_command
from phagemine.gui.services.inputs import (UploadedGenome, default_results_root, project_name,
                                            stage_uploads, validate_uploads)
from phagemine.gui.services.results import annotation_records, discovery_tables, downloadable_files, figures, result_zip
from phagemine.gui.services.runs import poll_run, start_run
from phagemine.gui.services.platform_support import windows_wsl_status
from phagemine.gui.services.status import (capability_rows, doctor_rows, platform_backend_note,
                                           read_run_status, workflow_readiness)
from phagemine.preflight import doctor
from phagemine.resources import EvidenceResourceManager, ResourceType

st.set_page_config(page_title="PhageMine", page_icon="🧬", layout="wide")
st.title("PhageMine")
st.caption(f"Local phage annotation and discovery · GUI v{GUI_VERSION} · Files never leave this computer")

PAGES = ("Home / New Analysis", "Environment / Database Status", "Run Monitor",
         "Annotation Results", "Discovery Results", "Figures", "Export / Downloads")
page = st.sidebar.radio("Workspace", PAGES, key="workspace")
result_dir = st.sidebar.text_input("Open result directory", st.session_state.get("result_dir", ""))
st.session_state.result_dir = result_dir
if os.environ.get("PHAGEMINE_DESKTOP") and st.sidebar.button("Quit PhageMine"):
    st.sidebar.success("PhageMine is closing. You may close this browser tab.")
    threading.Timer(0.5, lambda: os._exit(0)).start()


def caution() -> None:
    st.info("Unknown does not mean novel. NO_EXTERNAL_MATCH does not prove novelty. Genomic context supports but does not prove function. Computational predictions are not experimental confirmation.")


def current_doctor() -> dict:
    try:
        return doctor()
    except (OSError, ValueError, RuntimeError):
        return {"capabilities": {}, "executables": [], "resources": []}


if page == "Home / New Analysis":
    st.subheader("New analysis")
    st.markdown("### 1. Input")
    uploaded = st.file_uploader("Drag & drop FASTA files here or Browse files",
                                type=["fa", "fasta", "fna", "fas"], accept_multiple_files=True,
                                help="Choose one genome for a single analysis or two or more genomes for a cohort.")
    upload_items = [UploadedGenome(item.name, item.getvalue()) for item in uploaded]
    upload_rows, upload_error = [], None
    if upload_items:
        try:
            upload_rows = validate_uploads(upload_items)
            st.success(f"{len(upload_rows)} FASTA file{'s' if len(upload_rows) != 1 else ''} ready")
            st.dataframe(upload_rows, use_container_width=True, hide_index=True)
        except ValueError as exc:
            upload_error = str(exc); st.error(upload_error)

    advanced_source, source, forced_cohort = False, "", False
    with st.expander("Advanced input options"):
        advanced_source = st.checkbox("Use an existing local file or directory instead of uploads")
        if advanced_source:
            source_kind = st.radio("Local input type", ("Single FASTA file", "Directory of FASTA files"), horizontal=True)
            source = st.text_input("Local FASTA path" if source_kind.startswith("Single") else "Local cohort directory")
            forced_cohort = source_kind.startswith("Directory")

    cohort = forced_cohort if advanced_source else len(upload_items) > 1
    input_present = bool(source.strip()) if advanced_source else bool(upload_items)
    st.markdown("### 2. Analysis settings")
    col1, col2, col3 = st.columns(3)
    workflow_label = col1.selectbox("Workflow", ("Annotation", "Discovery", "Both"))
    mode = {"Annotation": "annotate", "Discovery": "discover", "Both": "both"}[workflow_label]
    evidence_options = ("core", "standard", "full") if cohort else ("core",)
    evidence = col2.selectbox("Evidence level", evidence_options,
                              help="Standard and full evidence profiles apply to cohort workflows and require validated databases.")
    threads = col3.number_input("Threads", 1, 256, min(4, os.cpu_count() or 1), 1)
    project = st.text_input("Project name", "phagemine-analysis")
    with st.expander("Advanced output location"):
        results_root = st.text_input("Results folder", str(default_results_root()))
    try:
        output = Path(results_root).expanduser() / project_name(project)
        output_error = None
    except ValueError as exc:
        output, output_error = Path(results_root).expanduser(), str(exc)
        st.error(output_error)
    st.caption(f"Results will be saved to: {output}")

    payload = current_doctor()
    ready, readiness_message = workflow_readiness(payload, mode=mode, evidence=evidence, cohort=cohort)
    (st.success if ready else st.warning)(readiness_message)
    resume = st.checkbox("Resume existing batch checkpoints", disabled=not cohort)
    preview_source = source if advanced_source else ("<uploaded-genome.fasta>" if len(upload_items) == 1 else "<uploaded-cohort-directory>")
    with st.expander("Advanced / command preview"):
        if input_present:
            try:
                preview = build_cli_args(preview_source, output, mode=mode, evidence=evidence,
                                         threads=int(threads), cohort=cohort, resume=resume)
                st.code(display_command(preview), language="shell")
            except ValueError as exc:
                st.error(str(exc))
        else:
            st.caption("Add FASTA input to preview the equivalent PhageMine command.")

    blocked = not input_present or bool(upload_error) or bool(output_error) or not ready
    if st.button("Start Analysis", type="primary", disabled=blocked):
        staged = None
        try:
            if advanced_source:
                actual_source = Path(source).expanduser()
                if not actual_source.exists():
                    raise ValueError("The selected local input does not exist")
                actual_cohort = forced_cohort
            else:
                staged = stage_uploads(upload_items)
                actual_source, actual_cohort = staged.path, staged.cohort
            actual = build_cli_args(actual_source, output, mode=mode, evidence=evidence,
                                    threads=int(threads), cohort=actual_cohort, resume=resume)
            handle = start_run(actual, input_cleanup=staged.root if staged else None)
            st.session_state.run_id = handle.run_id
            st.session_state.result_dir = str(output)
            st.session_state.workspace = "Run Monitor"
            st.rerun()
        except (OSError, ValueError) as exc:
            if staged:
                import shutil
                shutil.rmtree(staged.root, ignore_errors=True)
            st.error(f"Could not start analysis: {exc}")

elif page == "Environment / Database Status":
    st.subheader("Environment and database readiness")
    st.write(platform_backend_note())
    wsl = windows_wsl_status()
    if wsl["applicable"]:
        (st.success if wsl["state"] == "AVAILABLE" else st.warning)(f"WSL2 backend: {wsl['state']} — {wsl['detail']}")
    try:
        payload = doctor()
        st.markdown("### Analysis capabilities")
        st.dataframe(capability_rows(payload), use_container_width=True, hide_index=True)
        st.markdown("### Tools and evidence resources")
        st.dataframe(doctor_rows(payload), use_container_width=True, hide_index=True)
        not_ready = [row for row in doctor_rows(payload) if row["state"] not in {"READY", "OPTIONAL"}]
        if not_ready:
            st.warning("Some capabilities are unavailable. Install missing tools separately or register validated local evidence resources below.")
    except (OSError, ValueError, RuntimeError) as exc:
        st.error(f"Diagnostics failed: {exc}")

    with st.expander("Setup / Configure Evidence Databases"):
        st.write("Large databases are stored outside PhageMine so they can be updated independently. No database is downloaded automatically without a verified release URL and checksum.")
        kind = st.selectbox("Resource", ("PFAM", "VOGDB", "SWISSPROT", "PHROGS"))
        name = st.text_input("Registration name", kind)
        resource_path = st.text_input("Prepared database path")
        sidecar_label = {"VOGDB": "Annotation mapping path", "SWISSPROT": "Metadata path", "PHROGS": "Annotation mapping path"}.get(kind)
        sidecar = st.text_input(sidecar_label) if sidecar_label else ""
        if st.button("Register and validate resource"):
            try:
                provenance = {}
                if sidecar:
                    provenance["metadata_path" if kind == "SWISSPROT" else "annotations_path"] = sidecar
                manager = EvidenceResourceManager()
                manager.register(name, ResourceType(kind), resource_path, provenance=provenance)
                result = manager.validate(name)
                if result and result.get("status") == "READY": st.success(f"{name} is READY")
                else: st.error("Resource registered but not ready: " + "; ".join((result or {}).get("validation_errors", [])))
            except (OSError, ValueError) as exc:
                st.error(f"Could not register resource: {exc}")
    with st.expander("Setup / Fix scientific tools"):
        st.write("PhageMine validates tools already installed on this computer; it does not silently download or replace scientific executables.")
        st.markdown("- **macOS:** install PHANOTATE, HMMER, MMseqs2, and DIAMOND using trusted project packages, then restart PhageMine.\n- **Windows:** DIAMOND has a native binary and MMseqs2 has preview support, but upstream HMMER is POSIX-only. Use a maintained WSL2 scientific backend for full evidence workflows. Administrator approval and a restart may be required; PhageMine does not enable WSL without consent.")

elif page == "Run Monitor":
    st.subheader("Run monitor")
    run_id = st.session_state.get("run_id")
    live = poll_run(run_id) if run_id else None
    if live:
        st.metric("Process", live["state"])
        st.write("Output directory:", live["output_directory"])
        if live.get("command"):
            with st.expander("Executed command"): st.code(display_command(live["command"]), language="shell")
        if live.get("stderr"):
            with st.expander("Live diagnostic log", expanded=live["state"] == "FAILED"): st.text(live["stderr"][-20000:])
        if live["state"] == "FAILED": st.error(f"PhageMine exited with code {live['exit_code']}. No result was fabricated.")
        elif live["state"] == "COMPLETE": st.success("Analysis completed successfully.")
    status = read_run_status(result_dir) if result_dir else {"state": "NOT FOUND", "completed_stages": [], "warnings": [], "errors": []}
    if status.get("current_stage"): st.write("Current stage:", status["current_stage"])
    if status["completed_stages"]: st.write("Completed stages:", ", ".join(status["completed_stages"]))
    for warning in status["warnings"]: st.warning(warning)
    for error in status["errors"]: st.error(error)
    if status.get("samples"): st.dataframe(status["samples"], use_container_width=True)
    st.button("Refresh status")

elif page == "Annotation Results":
    st.subheader("Annotation results")
    try:
        records = annotation_records(result_dir) if result_dir else []
        if not records: st.warning("No native annotation outputs were found in this directory.")
        else:
            classification = st.multiselect("Functional classification", sorted({str(r.get("functional_state")) for r in records if r.get("functional_state")}))
            confidence = st.multiselect("Confidence", sorted({str(r.get("confidence") or r.get("functional_confidence")) for r in records if r.get("confidence") or r.get("functional_confidence")}))
            keyword = st.text_input("Keyword (protein ID or proposed function)").lower()
            filtered = [r for r in records if (not classification or r.get("functional_state") in classification)
                        and (not confidence or (r.get("confidence") or r.get("functional_confidence")) in confidence)
                        and (not keyword or keyword in f"{r.get('protein_id','')} {r.get('proposed_function','')} {r.get('annotation','')}".lower())]
            columns = ["protein_id", "sample", "start", "end", "strand", "functional_state", "proposed_function", "confidence", "functional_category", "pmf_id"]
            st.dataframe([{k: r.get(k) for k in columns} for r in filtered], use_container_width=True, hide_index=True)
            selected = st.selectbox("Protein detail", [r["protein_id"] for r in filtered]) if filtered else None
            if selected:
                record = next(r for r in filtered if r["protein_id"] == selected)
                st.json({k: v for k, v in record.items() if k not in {"sequence", "cds"}})
                st.text_area("Amino-acid sequence", record.get("sequence") or "Not available", height=130)
                st.text_area("CDS sequence", record.get("cds") or "Not available", height=130)
    except ValueError as exc: st.error(str(exc))

elif page == "Discovery Results":
    st.subheader("Discovery results"); caution()
    try:
        tables = discovery_tables(result_dir)
        table_name = st.selectbox("Native output table", tuple(tables))
        rows = tables[table_name]
        keyword = st.text_input("Filter rows")
        if keyword: rows = [r for r in rows if keyword.lower() in " ".join(map(str, r.values())).lower()]
        st.dataframe(rows, use_container_width=True, hide_index=True)
    except ValueError as exc: st.error(str(exc))

elif page == "Figures":
    st.subheader("Generated figures")
    found = figures(result_dir) if result_dir else []
    if not found: st.warning("No generated figures were found.")
    for path in found:
        st.caption(str(path.relative_to(Path(result_dir))))
        st.image(str(path), use_column_width=True)
        data_dir = path.parent.parent / "figure_data"
        if data_dir.is_dir(): st.write("Figure source data:", [p.name for p in sorted(data_dir.iterdir()) if p.is_file()])

else:
    st.subheader("Export and downloads")
    root = Path(result_dir)
    if not root.is_dir(): st.warning("The result directory does not exist.")
    else:
        files = downloadable_files(root)
        st.dataframe([{"file": str(p.relative_to(root)), "bytes": p.stat().st_size} for p in files], use_container_width=True, hide_index=True)
        st.download_button("Download complete native result directory as ZIP", result_zip(root), file_name=f"{root.name}.zip", mime="application/zip")
