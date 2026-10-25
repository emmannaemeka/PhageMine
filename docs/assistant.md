# Ask PhageMine: the conversational assistant

> Describe what you want to investigate. PhageMine will do the analysis and help you understand the results.

The **Ask PhageMine** page is the first page of `phagemine gui`. A researcher uploads FASTA files and types requests in ordinary English. A language model understands the request, chooses PhageMine tools, starts analyses, reads their outputs, and explains the actual results. PhageMine's existing, tested workflows do all scientific computation.

## Deploy in one day

**1. Install**

```bash
pip install "phagemine[assistant]"     # the extra is only needed for Anthropic models
phagemine databases install --all      # evidence databases; can take hours to download
phagemine doctor                       # confirm tools and databases
```

Install MMseqs2 (for example `conda install -c bioconda mmseqs2`). Without it, comparisons find homologs only through shared PHROG assignments and identical sequences; the assistant needs your explicit approval (an Approve button) before comparing that way. Install NCBI `table2asn` if you want official validation of submission packages.

**2. Choose a language model**

The assistant needs a language model to understand requests. Choose one:

| Option | Cost | Data leaves the computer? | Setup |
| --- | --- | --- | --- |
| **Ollama, local** (default when running) | Free | No | Install [Ollama](https://ollama.com), then `ollama pull qwen3:8b`. Keep Ollama open. |
| **Free cloud tier** (for example Google Gemini) | Free within daily limits | Yes; free tiers may use your data to improve their products | See below |
| **Anthropic Claude** | Paid per use; [research credits](https://support.claude.com/en/articles/11199177-anthropics-ai-for-science-program) available on application | Yes | `export ANTHROPIC_API_KEY=...` |

*Local model.* PhageMine finds a running Ollama automatically and picks the best installed model in this order: `qwen3:32b`, `qwen3:30b-a3b`, `qwen3:14b`, `qwen3:8b`, `qwen3:4b`. Bigger is more reliable at multi-step requests but needs more memory: `qwen3:8b` runs on a 16 GB laptop (slowly without a GPU); `qwen3:14b` or larger is better with a GPU. Override with `PHAGEMINE_ASSISTANT_MODEL=qwen3:14b`.

*Free cloud tier (OpenAI-compatible).* For Google Gemini, create a free key in Google AI Studio, then:

```bash
export PHAGEMINE_ASSISTANT_PROVIDER=openai-compatible
export PHAGEMINE_ASSISTANT_BASE_URL=https://generativelanguage.googleapis.com/v1beta/openai
export PHAGEMINE_ASSISTANT_MODEL=<a free-tier Flash model ID from AI Studio>
export PHAGEMINE_ASSISTANT_API_KEY=<your Gemini key>
```

Other services with this format work the same way (OpenRouter: `https://openrouter.ai/api/v1`; Groq: `https://api.groq.com/openai/v1`). One request to the assistant can make several model calls, so low daily limits run out quickly.

**3. Start**

```bash
phagemine gui
```

The browser opens on **Ask PhageMine** and shows which model is in use. Create a project in the sidebar, add FASTA files, and type a request.

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
| Tools | `assistant/tools.py` | 17 tools. Arguments are checked at run time against each schema (types, allowed values, bounds). Read tools answer from existing results; job tools start work and return a job ID only. No shell or code-execution tool exists. |
| Job worker | `assistant/jobs.py`, `assistant/runner.py` | Separate process per job. **COMPLETE requires exit code 0 and validated outputs**: annotations pass PhageMine's `validate_run` (protein IDs, sequences and coordinates agree across files); comparisons are parsed and every protein and source run is checked against its manifest. Otherwise FAILED with the reason. State changes are locked, and a final status is never overwritten. Dead workers are detected, including before they start. |
| Project store | `assistant/store.py` | One folder per project; inputs copied read-only; each job writes a new `runs/<job_id>/` folder. Tools cannot reach paths outside the project. |
| Agent loop | `assistant/agent.py` | Model, tools, model, until an answer. A compact project ledger is rebuilt every turn. |
| Grounding check | `assistant/grounding.py` | Before display, against tool results only (what the researcher typed is not evidence): protein and job IDs and numbers of 10 or more must appear in tool results; "completed" must match the status of the job or sample named in the same sentence; "novel" is refused. One automatic retry; if it fails again the written answer is withheld and the tool results are shown instead. |
| Approvals | `assistant/approvals.py` | Actions that need the researcher's decision (comparing without MMseqs2, applying NCBI corrections) create a request shown with Approve / Decline buttons. Only the buttons decide; tools and the model cannot. |
| NCBI support | `assistant/ncbi.py` | Error triage, guarded corrections, `phagemine revise`. |
| Model | `assistant/llm.py` | Ollama (local, free), any OpenAI-compatible service, or Anthropic, behind one interface. |
| Page | `gui/assistant_page.py` | Chat, uploads, job panel, downloads. |

### Comparisons

`start_comparison` runs `phagemine compare` once per genome as the query. Results are merged into deduplicated protein pairs, each marked with the direction(s) that found it, so counts and tables agree and input order does not change the answer. Results always come from the exact annotation runs the comparison used, even after re-annotation.

`phagemine compare` now writes `comparison_manifest.json` recording its source runs (with input checksums), thresholds, and how many MMseqs2 searches were attempted, succeeded and failed. The assistant reports MMseqs2 as used only when searches actually succeeded; otherwise it states that homologs came only from shared PHROG assignments and identical sequences, so counts are a lower bound. Comparing without MMseqs2 requires the researcher's Approve button.

### NCBI submission and error correction

1. The assistant reads the package PhageMine already wrote (`genbank_submission/`) and lists missing metadata.
2. Missing source details (host, collection date, location, isolation source, submitter, authors) are **asked for**. A proposed value must appear in what the researcher typed, but that alone is not confirmation: every value is shown again in the approval card. Metadata already in the original package is preserved.
3. Pasted errors are classified:
   - `FIXABLE_LABEL`: product, note, partial flag, locus tag;
   - `FIXABLE_METADATA`: source and submitter details;
   - `NEEDS_REANALYSIS`: internal stops, frameshifts, start/stop codons, overlaps, short CDSs, gaps. These need gene-model review; relabelling cannot fix them;
   - `UNKNOWN`: shown as received, no fix proposed.
4. The assistant shows a preview and the page shows an **Approve / Decline** card listing every correction the new package will contain. **Only the Approve button authorises it; a typed message never does.** The approval is bound to a checksum of those exact corrections and to the source annotation, so any later change voids it.
5. `phagemine revise` writes a new package with a correction audit and re-runs `table2asn` when installed. Corrections accumulate across rounds.
6. Product names may only become "hypothetical protein", be re-formatted, or gain an uncertainty qualifier such as "putative". Removing a qualifier or naming a function more specifically needs evidence review and curation; words inside evidence descriptions never authorise a name. Coordinates and partial flags are curation decisions and cannot be changed this way.
7. Corrections from successive rounds are merged field by field, so a later note does not erase an earlier product fix. A finished revision job means a new package was written and validated; it does not mean NCBI will accept it.

PhageMine never submits anything to NCBI.

## Data sent to the model provider

With a local Ollama model nothing leaves the computer. With a cloud model, each turn sends your messages, a project summary, and tool results: file names, gene IDs, coordinates, product names, database hits and scores, and NCBI messages. Genome and protein sequences are not sent unless you ask to see one. Every model call's complete payload (system prompt, messages and tool list, exactly as sent) is logged in the project's `outbound.jsonl`. Check your provider's current data-retention terms before using unpublished genomes.

## Limitations

- Tested with a scripted stand-in model and in demo mode. Measure tool-selection accuracy with your chosen model on real genomes before relying on it.
- The grounding check verifies IDs, numbers of 10 or more, completion claims and novelty wording. It does not verify the meaning of every sentence or bind each claim to its evidence field; evidence tables remain the authority. Small numbers (0-9) are not checked, and the novelty rule is a word rule, not a contextual judgement.
- Not yet built (from external review): comparison semantics owned by the core (canonical pairs, reused search databases), claim-level binding of every statement to a result field, multi-user server mode with authentication and quotas, small-model mode, and a pre-registered evaluation of about 200 scenarios with real models.
- The NCBI error catalogue is keyword-based. Extend `CATALOGUE` in `assistant/ncbi.py` with messages from your real submissions.
- One researcher per project. Concurrent annotations share the machine's CPUs.
- Local models are less reliable than large cloud models at multi-step requests. Measure with your model before relying on it.

## Tests

`tests/test_assistant_review.py` holds 22 counterexamples from external review, each asserting the corrected behaviour. `tests/test_assistant.py` covers job success and failure rules, input protection, project isolation, both-direction comparison, the figure's source table, NCBI guards, the grounding check and the agent loop. `tests/test_gui_app.py::test_conversational_analysis_end_to_end` drives the real page through annotation, comparison, the automatic report, and a follow-up question.
