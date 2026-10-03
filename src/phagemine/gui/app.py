"""PhageMine GUI v0.1 — a local Streamlit view over native CLI artifacts."""
from __future__ import annotations

import io
import tempfile
import zipfile
from pathlib import Path

import streamlit as st

from phagemine.gui import GUI_VERSION
from phagemine.gui.services.execution import build_cli_args, display_command, execute, save_uploads, validate_input
from phagemine.gui.services.results import annotation_records, discovery_tables, downloadable_files, export_files, figures
from phagemine.gui.services.status import doctor_rows, read_run_status
from phagemine.preflight import doctor

st.set_page_config(page_title="PhageMine GUI", page_icon="🧬", layout="wide")
st.title("PhageMine GUI")
st.caption(f"Local scientific interface · GUI v{GUI_VERSION} · Same PhageMine analysis engine")

PAGES = ("Home / New Analysis", "Environment / Database Status", "Run Monitor",
         "Annotation Results", "Discovery Results", "Figures", "Export / Downloads")
page = st.sidebar.radio("Workspace", PAGES)
result_dir = st.sidebar.text_input("Result directory", key="result_dir", value="results/gui_run")


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
    if not cohort:
        st.caption("The single-genome CLI uses the available registered evidence resources. Explicit standard/full profiles are available for cohort inputs.")
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
        else:
            try:
                with tempfile.TemporaryDirectory(prefix="phagemine-gui-") as temporary:
                    actual_source = source
                    if uploads:
                        actual_source = str(save_uploads(uploads, Path(temporary)))
                    actual_source = str(validate_input(actual_source, cohort=cohort))
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
                            # A widget's state cannot be changed after it was rendered.
                            st.success(f"Results saved to {Path(output).expanduser()}. Open this path in the sidebar to review them.")
            except (OSError, ValueError) as exc:
                st.error(f"Could not start analysis: {exc}")

elif page == "Environment / Database Status":
    st.subheader("Environment and database status")
    st.caption("States come from PhageMine doctor and resource validation; paths alone are not treated as readiness.")
    deep = st.checkbox("Check operational database formats", value=False)
    try:
        rows = doctor_rows(doctor(deep=deep))
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
            selected = st.selectbox("Protein detail", [r["record_id"] for r in filtered]) if filtered else None
            if selected:
                record = next(r for r in filtered if r["record_id"] == selected)
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
    root = Path(result_dir).expanduser()
    if not root.is_dir(): st.warning("The result directory does not exist.")
    else:
        files = downloadable_files(root)
        st.dataframe([{"file": str(p.relative_to(root)), "bytes": p.stat().st_size} for p in files], use_container_width=True, hide_index=True)
        if st.button("Prepare result ZIP"):
            buffer = io.BytesIO()
            resolved_root = root.resolve()
            with zipfile.ZipFile(buffer, "w", zipfile.ZIP_DEFLATED) as archive:
                for path in export_files(root): archive.write(path, path.relative_to(resolved_root))
            st.download_button("Download result directory as ZIP", buffer.getvalue(), file_name=f"{root.name}.zip", mime="application/zip")
