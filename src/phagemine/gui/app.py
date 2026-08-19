"""PhageMine GUI v0.1 — a local Streamlit view over native CLI artifacts."""
from __future__ import annotations

import io
import tempfile
import zipfile
from pathlib import Path

import streamlit as st

from phagemine.gui import GUI_VERSION
from phagemine.gui.services.execution import build_cli_args, display_command, execute
from phagemine.gui.services.results import annotation_records, discovery_tables, downloadable_files, figures
from phagemine.gui.services.status import doctor_rows, read_run_status
from phagemine.preflight import doctor

st.set_page_config(page_title="PhageMine GUI", page_icon="🧬", layout="wide")
st.title("PhageMine GUI")
st.caption(f"Local scientific interface · GUI v{GUI_VERSION} · Same PhageMine analysis engine")

PAGES = ("Home / New Analysis", "Environment / Database Status", "Run Monitor",
         "Annotation Results", "Discovery Results", "Figures", "Export / Downloads")
page = st.sidebar.radio("Workspace", PAGES)
result_dir = st.sidebar.text_input("Result directory", st.session_state.get("result_dir", "results/gui_run"))
st.session_state.result_dir = result_dir


def caution() -> None:
    st.info("Unknown does not mean novel. NO_EXTERNAL_MATCH does not prove novelty. Genomic context supports but does not prove function. Computational predictions are not experimental confirmation.")


if page == "Home / New Analysis":
    st.subheader("New analysis")
    source_kind = st.radio("Input", ("Single FASTA path", "Cohort directory", "Upload FASTA files"), horizontal=True)
    source = ""
    uploads = None
    if source_kind == "Single FASTA path":
        source = st.text_input("FASTA genome path")
    elif source_kind == "Cohort directory":
        source = st.text_input("Directory containing FASTA genomes")
    else:
        uploads = st.file_uploader("FASTA genome files", type=["fa", "fasta", "fna", "fas"], accept_multiple_files=True)
    col1, col2, col3 = st.columns(3)
    mode = col1.selectbox("Workflow", ("Annotate", "Discover", "Both")).lower()
    cohort = source_kind != "Single FASTA path"
    evidence_choices = ("core", "standard", "full") if cohort else ("core",)
    evidence = col2.selectbox("Evidence level", evidence_choices)
    threads = col3.number_input("Threads", 1, 256, min(4, 256), 1)
    output = st.text_input("Output directory", result_dir)
    resume = st.checkbox("Resume existing batch checkpoints", disabled=not cohort)

    preview_source = source or ("<uploaded-fasta-directory>" if uploads else "<input>")
    try:
        args = build_cli_args(preview_source, output, mode=mode, evidence=evidence,
                              threads=int(threads), cohort=cohort, resume=resume)
        st.code(display_command(args), language="shell")
    except ValueError as exc:
        st.error(str(exc)); args = []
    if st.button("Start analysis", type="primary"):
        if not output.strip():
            st.error("Choose a valid output directory.")
        elif source_kind == "Upload FASTA files" and not uploads:
            st.error("Upload at least one FASTA file.")
        elif source_kind != "Upload FASTA files" and not Path(source).expanduser().exists():
            st.error("The selected input does not exist.")
        else:
            try:
                with tempfile.TemporaryDirectory(prefix="phagemine-gui-") as temporary:
                    actual_source = source
                    if uploads:
                        upload_dir = Path(temporary)
                        for upload in uploads:
                            (upload_dir / Path(upload.name).name).write_bytes(upload.getvalue())
                        actual_source = str(upload_dir)
                    actual = build_cli_args(actual_source, output, mode=mode, evidence=evidence,
                                            threads=int(threads), cohort=cohort, resume=resume)
                    with st.status("PhageMine analysis running", expanded=True) as box:
                        st.code(display_command(actual), language="shell")
                        result = execute(actual)
                        if result.stdout: st.text(result.stdout)
                        if result.stderr: st.text(result.stderr)
                        if result.exit_code:
                            box.update(label="Analysis failed", state="error")
                            st.error(f"PhageMine exited with code {result.exit_code}. See diagnostic output above.")
                        else:
                            box.update(label="Analysis complete", state="complete")
                            st.session_state.result_dir = output
            except (OSError, ValueError) as exc:
                st.error(f"Could not start analysis: {exc}")

elif page == "Environment / Database Status":
    st.subheader("Environment and database status")
    st.caption("States come from PhageMine doctor and resource validation; paths alone are not treated as readiness.")
    try:
        rows = doctor_rows(doctor())
        st.dataframe(rows, use_container_width=True, hide_index=True)
    except (OSError, ValueError, RuntimeError) as exc:
        st.error(f"Diagnostics failed: {exc}")

elif page == "Run Monitor":
    st.subheader("Run monitor")
    status = read_run_status(result_dir)
    st.metric("State", status["state"])
    st.write("Output directory:", status["output_directory"])
    if status.get("current_stage"): st.write("Current stage:", status["current_stage"])
    if status["completed_stages"]: st.write("Completed stages:", ", ".join(status["completed_stages"]))
    for warning in status["warnings"]: st.warning(warning)
    for error in status["errors"]: st.error(error)
    if status.get("samples"): st.dataframe(status["samples"], use_container_width=True)
    if status.get("resume_available"):
        st.caption("Resume is available from Home for batch runs and delegates to PhageMine's existing --resume-existing checkpoints.")
    st.button("Refresh status")

elif page == "Annotation Results":
    st.subheader("Annotation results")
    try:
        records = annotation_records(result_dir)
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
    found = figures(result_dir)
    if not found: st.warning("No generated figures were found.")
    for path in found:
        st.caption(str(path.relative_to(Path(result_dir))))
        st.image(str(path), use_column_width=True)
        data_dir = path.parent.parent / "figure_data"
        if data_dir.is_dir():
            st.write("Figure source data:", [p.name for p in sorted(data_dir.iterdir()) if p.is_file()])

else:
    st.subheader("Export and downloads")
    root = Path(result_dir)
    if not root.is_dir(): st.warning("The result directory does not exist.")
    else:
        files = downloadable_files(root)
        st.dataframe([{"file": str(p.relative_to(root)), "bytes": p.stat().st_size} for p in files], use_container_width=True, hide_index=True)
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
            for path in sorted(p for p in root.rglob("*") if p.is_file()): archive.write(path, path.relative_to(root))
        st.download_button("Download completed result directory as ZIP", buffer.getvalue(), file_name=f"{root.name}.zip", mime="application/zip")
