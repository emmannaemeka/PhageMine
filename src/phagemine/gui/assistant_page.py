"""Ask PhageMine: the conversational assistant page.

All logic lives in ``phagemine.assistant``; this page only displays the
conversation, uploads files, and polls background jobs.
"""
from __future__ import annotations

import json

import streamlit as st

from phagemine.assistant import agent, approvals, jobs
from phagemine.assistant.llm import ProviderUnavailable, default_provider
from phagemine.assistant.store import Project, projects_root
from phagemine.assistant.tools import demo_mode

HEADLINE = "Describe what you want to investigate. PhageMine will do the analysis and help you understand the results."
EXAMPLES = ("Here are two phage genomes. Annotate them, compare them, and show me which genes they share.",
            "Which of these genes have unknown functions?",
            "What evidence supports the annotation of PM_000012?",
            "Make a figure showing the comparison.",
            "Prepare this genome for NCBI submission.")


def _project() -> Project | None:
    projects = Project.list()
    st.sidebar.markdown("**Assistant project**")
    names = {f"{p.name} ({p.id})": p for p in projects}
    pending = st.session_state.pop("assistant_project_pending", None)
    if pending in names:  # a widget's value can only be set before the widget is drawn
        st.session_state.assistant_project = pending
    choice =st.sidebar.selectbox("Project", list(names) or ["No projects yet"], key="assistant_project",
                                  help="Each project keeps its own files, analyses and conversation.")
    with st.sidebar.form("new_project", clear_on_submit=True):
        name = st.text_input("New project name")
        if st.form_submit_button("Create project") and name.strip():
            created = Project.create(name.strip())
            st.session_state.assistant_project_pending = f"{created.name} ({created.id})"
            st.rerun()
    st.sidebar.caption(f"Projects are stored in {projects_root()}")
    return names.get(choice)


def _provider():
    if "assistant_provider" not in st.session_state:
        try:
            st.session_state.assistant_provider = default_provider()
            st.session_state.assistant_provider_error = None
        except ProviderUnavailable as exc:
            st.session_state.assistant_provider = None
            st.session_state.assistant_provider_error = str(exc)
    return st.session_state.assistant_provider


def _files(project: Project, names: list[str]) -> None:
    for name in names:
        path = project.root / ("figures" if (project.root / "figures" / name).exists() else "summaries") / name
        if not path.is_file():
            continue
        if path.suffix == ".png":
            st.image(str(path), caption=name)
        st.download_button(f"Download {name}", path.read_bytes(), file_name=name, key=f"dl-{project.id}-{name}-{id(names)}")


def _transcript(project: Project) -> None:
    for item in project.transcript():
        if item["role"] == "user":
            with st.chat_message("user"):
                st.markdown(item["text"])
        elif item["role"] == "event":
            st.caption("PhageMine: background analysis finished; the assistant reported on it below.")
        else:
            with st.chat_message("assistant"):
                st.markdown(item["text"] or "_(no answer)_")
                calls = item.get("tool_calls") or []
                if calls:
                    with st.expander(f"What PhageMine did ({len(calls)} steps)"):
                        for call in calls:
                            mark = "done" if call.get("ok") else f"error: {call.get('error')}"
                            st.markdown(f"- `{call['name']}` {json.dumps(call.get('input') or {})} — {mark}")
                _files(project, item.get("files") or [])


def _decide(project: Project, approval_id: str, approve: bool) -> None:
    """Button callback: the only code path that grants or declines an approval."""
    record = (approvals.approve if approve else approvals.decline)(project, approval_id)
    word = "approved" if approve else "declined"
    preview = record["payload"].get("id") if isinstance(record.get("payload"), dict) else None
    st.session_state.assistant_event = (f"[The researcher pressed the button and {word}: {record['title']} "
                                        f"(approval {record['id']}{', preview ' + preview if preview else ''}). Continue accordingly.]")


def _approvals(project: Project) -> None:
    for record in approvals.pending(project):
        with st.container(border=True):
            st.markdown(f"**Your approval is needed:** {record['title']}")
            payload = record.get("payload") or {}
            corrections = payload.get("corrections") if isinstance(payload, dict) else None
            if corrections:
                rows = [{"what": f"metadata: {k}", "new value": str(v)} for k, v in corrections.get("metadata", {}).items()]
                rows += [{"what": f"{pid}: {field}", "new value": str(value)}
                         for pid, fields in corrections.get("annotations", {}).items()
                         for field, value in fields.items() if field not in {"reason", "source"}]
                st.table(rows)
                st.caption("These are all the corrections the new package will contain, including earlier approved ones.")
            cols = st.columns(2)
            cols[0].button("Approve", key=f"approve-{record['id']}", type="primary",
                           on_click=_decide, args=(project, record["id"], True))
            cols[1].button("Decline", key=f"decline-{record['id']}", on_click=_decide, args=(project, record["id"], False))


@st.fragment(run_every=3)
def _jobs_panel(project: Project) -> None:
    records = [jobs.status(project, j["id"]) for j in project.jobs()]
    active = [r for r in records if r["status"] not in jobs.FINAL]
    if active:
        st.markdown("**Running analyses**")
        for record in active:
            label = record.get("sample") or ", ".join(record.get("samples", []))
            progress = record.get("progress")
            detail = f"{progress['stage']} ({progress['completed']}/{progress['total']})" if progress else record["status"].lower()
            cols = st.columns([5, 1])
            cols[0].write(f"{record['kind'].title()} · {label} · {detail}")
            if cols[1].button("Cancel", key=f"cancel-{record['id']}"):
                jobs.cancel(project, record["id"])
            if progress:
                cols[0].progress(min(progress["completed"] / max(progress["total"], 1), 1.0))
    if not st.session_state.get("assistant_event"):
        text = agent.job_event(project)
        if text:
            st.session_state.assistant_event = text
            st.rerun(scope="app")


def render() -> None:
    st.subheader("Ask PhageMine")
    st.write(HEADLINE)
    project = _project()
    provider = _provider()
    if demo_mode():
        st.warning("Demo mode: analyses use PhageMine's fixture gene caller and mock evidence. Results are not biology.")
    if provider is None:
        st.error(st.session_state.assistant_provider_error)
    else:
        local = str(getattr(provider, "name", "")).startswith("ollama:")
        st.caption(f"Model: {provider.name}" + (" · runs on this computer; nothing is sent online" if local else ""))
    if project is None:
        st.info("Create a project in the sidebar to begin.")
        return
    with st.expander("Data and privacy", expanded=False):
        st.markdown("Analyses run on this computer. With a local Ollama model, nothing leaves it. With a cloud model, "
                    "to understand your requests the assistant sends your messages "
                    "and tool results (file names, gene IDs, coordinates, product names, database hits and scores, "
                    "NCBI messages) to the language model provider. Genome and protein sequences are not sent unless "
                    "you ask to see one. Every payload is logged in `outbound.jsonl` in the project folder. "
                    "Your uploaded files are copied into the project and never modified.")
    uploads = st.file_uploader("Add genome files (FASTA)", accept_multiple_files=True,
                               type=["fa", "fasta", "fna", "fas", "fsa"], key=f"upload-{project.id}")
    for upload in uploads or []:
        marker = f"added-{project.id}-{upload.name}-{upload.size}"
        if marker not in st.session_state:
            try:
                project.add_input(upload.name, upload.getvalue())
                st.session_state[marker] = True
            except ValueError as exc:
                st.error(str(exc))
    inputs = project.inputs()
    if inputs:
        st.caption("Files in this project: " + ", ".join(i["file"] for i in inputs))
    if not project.transcript():
        with st.expander("Examples of what you can ask", expanded=True):
            st.markdown("\n".join(f"- {e}" for e in EXAMPLES))
    _transcript(project)
    _approvals(project)
    _jobs_panel(project)
    outputs = project.outputs()
    if outputs:
        with st.sidebar.expander(f"Saved outputs ({len(outputs)})"):
            for path in outputs:
                st.download_button(path.name, path.read_bytes(), file_name=path.name, key=f"side-{project.id}-{path.name}")
    pending = st.session_state.pop("assistant_event", None)
    if pending and provider is not None:
        with st.chat_message("assistant"):
            with st.spinner("Reviewing finished analyses..."):
                agent.run_turn(project, pending, provider, event=True)
        st.rerun()
    prompt = st.chat_input("Describe what you want to investigate", disabled=provider is None, max_chars=8000)
    if prompt:
        with st.chat_message("user"):
            st.markdown(prompt)
        with st.chat_message("assistant"):
            with st.status("Working...", expanded=True) as box:
                def step(item: dict) -> None:
                    if item["type"] == "tool_start":
                        box.write(f"Running `{item['name']}`")
                    elif item["type"] == "tool_end" and not item["ok"]:
                        box.write(f"`{item['name']}` reported: {item['error']}")
                try:
                    agent.run_turn(project, prompt, provider, on_step=step)
                    box.update(label="Done", state="complete", expanded=False)
                except Exception as exc:  # network or provider errors must not lose the conversation
                    box.update(label="The assistant could not respond", state="error")
                    st.error(f"{type(exc).__name__}: {exc}")
                    st.stop()
        st.rerun()
