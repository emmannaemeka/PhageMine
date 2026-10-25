"""Plain-English entry point using the existing local workflow and result readers."""
from __future__ import annotations

import json
from pathlib import Path

import streamlit as st

from phagemine.gui.services.assistant import find_records, interpret, summarize
from phagemine.gui.services.results import annotation_records, discovery_tables


def _navigate(page: str) -> None:
    st.session_state.workspace_page = page


def render(result_dir: str) -> None:
    st.subheader("Ask PhageMine")
    st.write("Describe your question. Explore the evidence behind your annotations, or start a new analysis.")
    st.caption("First version · Offline request matching, not a generative AI model. Your files stay local. "
               "Search reads existing annotations; it does not perform a new sequence search.")
    with st.expander("Things you can ask", expanded=True):
        st.markdown('- “Annotate my genome”\n- “Find the capsid gene”\n- “Explain PM_000001”\n'
                    '- “Summarize my results”\n- “Compare my genomes”\n- “Check my installation”')
    query = st.chat_input("Ask about your genomes or results", max_chars=1000)
    scope = str(Path(result_dir).expanduser().resolve())
    if st.session_state.get("assistant_scope") != scope:
        st.session_state.assistant_scope = scope
        st.session_state.assistant_question = ""
    if query:
        st.session_state.assistant_question = query
    question = st.session_state.get("assistant_question", "")
    if not question:
        st.info("For existing analyses, select the Result directory in the sidebar. For a new genome, ask to annotate it.")
        return
    with st.chat_message("user"):
        st.write(question)
    request = interpret(question)
    with st.chat_message("assistant"):
        if request.action == "annotate":
            st.write("I can prepare an annotation using PhageMine’s existing workflow. "
                     "Upload your FASTA file on the next screen, choose a new output folder, and select Start analysis. "
                     "Available evidence depends on the databases installed on this computer.")
            st.button("Upload genome and start annotation", on_click=_navigate, args=("Home / New Analysis",))
        elif request.action == "environment":
            st.write("The environment check lists the tools and databases available for analysis.")
            st.button("Check tools and databases", on_click=_navigate, args=("Environment / Database Status",))
        elif request.action in {"search", "summary"}:
            try:
                records = annotation_records(result_dir)
                if request.action == "summary":
                    st.write(summarize(records))
                else:
                    matches = find_records(records, request.query)
                    st.write(f'Searching saved gene labels and annotations for: “{request.query}”.')
                    if not records:
                        st.warning("No annotation outputs found. Choose a completed result directory first.")
                    elif not matches:
                        st.info("No matching annotation was found. This does not establish that the gene is absent. "
                                "Try its protein ID, gene symbol, or another product name.")
                    else:
                        st.write(f"Found {len(matches)} matching annotation record(s). Predicted functions require scientific review.")
                        columns = ("record_id", "start", "end", "strand", "gene", "product", "proposed_function", "confidence", "evidence_sources", "best_evidence")
                        st.dataframe([{k: row.get(k) for k in columns} for row in matches], hide_index=True)
                        selected = st.selectbox("Inspect supporting evidence", [r["record_id"] for r in matches])
                        row = next(r for r in matches if r["record_id"] == selected)
                        st.caption("Source directory: " + row["run_directory"])
                        st.json({k: v for k, v in row.items() if k not in {"sequence", "cds"}})
                        st.download_button("Download matching annotations", json.dumps(matches, indent=2), "phagemine-matches.json", "application/json")
                st.caption("Read from: " + scope + ". Values are taken directly from saved PhageMine outputs.")
            except (OSError, ValueError) as exc:
                st.error(f"Could not read the selected outputs: {exc}")
        elif request.action == "compare":
            try:
                tables = discovery_tables(result_dir)
                if tables["families"]:
                    st.write("Here are the protein-family results already computed for this output folder. "
                             "Use membership and recurrence tables to inspect which samples share a family. "
                             "I have not run a new comparison.")
                    for name in ("families", "members", "recurrence"):
                        with st.expander(name.title(), expanded=name == "families"):
                            st.dataframe(tables[name], hide_index=True)
                else:
                    st.info("No computed family comparison was found here. Use New Analysis with a cohort of genomes "
                            "and the Both workflow, then select the completed result directory.")
                st.button("Open comparison workflow", on_click=_navigate, args=("Home / New Analysis",))
            except (OSError, ValueError) as exc:
                st.error(str(exc))
        else:
            st.write("I could not confidently map that request to a supported operation. "
                     "Try one question at a time using an example above. For a gene lookup, include its name or protein ID. "
                     "This first version cannot answer unrestricted biological questions.")
