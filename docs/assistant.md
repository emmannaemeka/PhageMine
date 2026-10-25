# Ask PhageMine: the conversational assistant

> Describe what you want to investigate. PhageMine will do the analysis and help you understand the results.

The **Ask PhageMine** page is the first page of `phagemine gui`. A researcher uploads FASTA files and types requests in ordinary English. A language model understands the request, chooses PhageMine tools, starts analyses, reads their outputs, and explains the actual results. PhageMine's existing, tested workflows do all scientific computation.

## Deploy in one day

**1. Install (10 minutes)**

```bash
pip install "phagemine[assistant]"
phagemine databases install --all      # evidence databases; can take hours to download
phagemine doctor                       # confirm tools and databases
```

Install MMseqs2 (for example `conda install -c bioconda mmseqs2`). Without it, comparisons can only detect identical proteins, and the assistant will ask before proceeding that way. Install NCBI `table2asn` if you want official validation of submission packages.

**2. Set the model key**

```bash
export ANTHROPIC_API_KEY=...                      # from console.anthropic.com
export PHAGEMINE_ASSISTANT_MODEL=claude-sonnet-5-5   # optional; this is the default
export PHAGEMINE_PROJECTS=~/PhageMine/projects       # optional; this is the default
```

**3. Start**

```bash
phagemine gui
```

The browser opens on **Ask PhageMine**. Create a project in the sidebar, add FASTA files, and type a request.

**4. Try it safely first (optional)**

```bash
PHAGEMINE_ASSISTANT_DEMO=1 phagemine gui
```

Demo mode uses PhageMine's fixture gene caller and mock evidence, so a full conversation runs in seconds without databases. Every answer and saved summary is labelled as demonstration data.

**Lab server.** The GUI binds to `127.0.0.1` on purpose: it can read files and launch analyses. To use a shared workstation, keep it local there and connect with an SSH tunnel: `ssh -L 8501:127.0.0.1:8501 user@server`, then open `http://127.0.0.1:8501`. Do not expose it publicly without an authentication layer.

## What researchers can ask

- "Here are two phage genomes. Annotate them, compare them, and show me which genes they share."
- "Which of these genes have unknown functions?"
- "What evidence supports the annotation of PM_000014 in phageB?"
- "Make a figure showing the comparison."
- "Explain this more simply and save a summary."
- "Prepare this genome for NCBI submission." Then paste the validator or GenBank staff messages: "NCBI sent these errors: ..."

Analyses run in the background. The page shows progress, and when they end the assistant reports the results without being asked again. Closing the browser does not stop a running analysis.

## How it is built

| Part | Module | Role |
| --- | --- | --- |
| Tools | `assistant/tools.py` | 17 typed tools. Read tools answer from existing results; job tools start work and return a job ID only. No shell or code-execution tool exists. |
| Job worker | `assistant/jobs.py`, `assistant/runner.py` | Separate process per job. **COMPLETE requires exit code 0 and validated outputs**; otherwise FAILED with the reason. |
| Project store | `assistant/store.py` | One folder per project; inputs copied read-only; each job writes a new `runs/<job_id>/` folder. Tools cannot reach paths outside the project. |
| Agent loop | `assistant/agent.py` | Model, tools, model, until an answer. A compact project ledger is rebuilt every turn. |
| Grounding check | `assistant/grounding.py` | Before display: protein and job IDs must appear in tool results, "completed" requires a COMPLETE job, and "novel" is refused. One automatic retry, then a visible warning. |
| NCBI support | `assistant/ncbi.py` | Error triage, guarded corrections, `phagemine revise`. |
| Model | `assistant/llm.py` | Anthropic provider; any provider returning the same message shape can replace it. |
| Page | `gui/assistant_page.py` | Chat, uploads, job panel, downloads. |

### Comparisons

`start_comparison` runs `phagemine compare` once per genome as the query, so shared genes are counted from both sides and input order does not change the answer. The method used (MMseqs2 or exact sequence identity) is reported with every result.

### NCBI submission and error correction

1. The assistant reads the package PhageMine already wrote (`genbank_submission/`) and lists missing metadata.
2. Missing source details (host, collection date, location, isolation source, submitter, authors) are **asked for**. A metadata value is accepted only if the researcher typed it in the conversation.
3. Pasted errors are classified:
   - `FIXABLE_LABEL`: product, note, partial flag, locus tag;
   - `FIXABLE_METADATA`: source and submitter details;
   - `NEEDS_REANALYSIS`: internal stops, frameshifts, start/stop codons, overlaps, short CDSs, gaps. These need gene-model review; relabelling cannot fix them;
   - `UNKNOWN`: shown as received, no fix proposed.
4. The assistant shows a before/after preview. **Corrections are applied only after the researcher approves in a later message.**
5. `phagemine revise` writes a new package with a correction audit and re-runs `table2asn` when installed. Corrections accumulate across rounds.
6. Product names may become less specific (for example "hypothetical protein") but never more specific than the protein's recorded evidence. Coordinates cannot be changed this way.

PhageMine never submits anything to NCBI.

## Data sent to the model provider

Each turn sends your messages, a project summary, and tool results: file names, gene IDs, coordinates, product names, database hits and scores, and NCBI messages. Genome and protein sequences are not sent unless you ask to see one. Every outbound payload is logged in the project's `outbound.jsonl`. Check your provider's current data-retention terms before using unpublished genomes.

## Limitations

- Tested with a scripted stand-in model and in demo mode. Measure tool-selection accuracy with your chosen model on real genomes before relying on it.
- The grounding check verifies IDs, completion claims and novelty wording. It does not verify every number or sentence; evidence tables remain the authority.
- The NCBI error catalogue is keyword-based. Extend `CATALOGUE` in `assistant/ncbi.py` with messages from your real submissions.
- One researcher per project. Concurrent annotations share the machine's CPUs.
- No local-model provider is included yet.

## Tests

`tests/test_assistant.py` covers job success and failure rules, input protection, project isolation, both-direction comparison, the figure's source table, NCBI guards, the grounding check and the agent loop. `tests/test_gui_app.py::test_conversational_analysis_end_to_end` drives the real page through annotation, comparison, the automatic report, and a follow-up question.
