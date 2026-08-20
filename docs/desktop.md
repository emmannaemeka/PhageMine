# PhageMine Desktop release-candidate guide

PhageMine Desktop packages the existing Streamlit interface and PhageMine Python
engine with PyInstaller. It is a local browser-backed desktop application: the
windowless launcher binds Streamlit to `127.0.0.1`, chooses an available local
port, and opens the default browser using Python's cross-platform browser API.
It never binds to all interfaces and does not use a cloud service.

The desktop package does not implement a second annotation pipeline. Analysis
still runs through the same `phagemine run`, `phagemine mine`, and
`phagemine batch` command semantics used by the CLI.

The release-candidate package version is `1.1.0rc1`; its compatible future tag
name is `v1.1.0-rc.1`. This tag has not been created or published.

## PhageMine Desktop — macOS

1. Download the DMG matching the Mac: `PhageMine-macOS-AppleSilicon-1.1.0rc1.dmg`
   for arm64 Macs or `PhageMine-macOS-Intel-1.1.0rc1.dmg` for Intel x86_64 Macs.
2. Open the DMG and move `PhageMine.app` to Applications.
3. Launch **PhageMine** from Applications.
4. Review Environment / Database Status before the first analysis.

Release-candidate DMGs are not yet signed or notarized. They are CI validation
artifacts, not a published v1.1.0 release.

The two DMGs are separate native builds, not a claimed universal application.
The Intel build is produced on GitHub's native `macos-15-intel` runner with a
minimum supported deployment target of macOS 12.0. CI audits every embedded
Mach-O executable, dynamic library, and extension module for x86_64 architecture
and a declared minimum macOS version no later than 12.0 before upload.

## PhageMine Desktop — Windows

1. Download the `PhageMine-Windows-<commit>` artifact from the native Windows build.
2. Run `PhageMine-Windows-Setup.exe`.
3. Launch **PhageMine** from the Start menu or optional desktop shortcut.

The executable is built with the PyInstaller windowed mode, so no console stays
open during normal use. Release-candidate installers are not yet code-signed.

### Windows scientific-backend reality

The Windows desktop is a **technical preview**. The installer and frozen CLI/UI
are tested on a clean native Windows CI runner, but production core annotation
is deliberately fail-closed there:

| Tool | macOS | Native Windows | Distribution decision |
|---|---|---|---|
| PHANOTATE 1.6.7 | Bundled separate backend; GPL-3.0-or-later | Not bundled: upstream `fastpath` has no Windows wheel and includes POSIX `getopt.h` | macOS includes attribution, licenses, and exact Corresponding Source archives; Windows reports NOT READY |
| HMMER / hmmscan | Upstream POSIX support | Upstream does not support native Windows | Use a maintained WSL2 backend for full evidence; never substitute an algorithm |
| MMseqs2 | Official universal build | Official Windows preview; upstream recommends WSL2 | External validated install |
| DIAMOND | Official binary | Official binary; Visual C++ runtime required | External validated install |
| table2asn | Optional | Optional | External optional install |

Authoritative project references:

- PHANOTATE: <https://github.com/deprekate/PHANOTATE>
- HMMER: <https://hmmer.org/>
- MMseqs2: <https://github.com/soedinglab/MMseqs2>
- DIAMOND installation: <https://github.com/bbuchfink/diamond/wiki/2.-Installation>
- NCBI table2asn: <https://www.ncbi.nlm.nih.gov/genbank/table2asn/>

PhageMine displays WSL2 availability separately from `phagemine doctor` output.
It does not install or enable WSL2 automatically because that can require
administrator approval and a restart. Until a Windows backend is validated by
doctor, the GUI must show **NOT READY** and block affected workflows.

The macOS PHANOTATE backend is an unmodified upstream scientific program run
as a separate process. The macOS application carries PHANOTATE, fastpath, and genbank
license texts plus their exact pinned source archives. No gene-calling
algorithm, default, or threshold is substituted. Other scientific executables
and all large evidence databases remain external and are never marked ready
unless the existing doctor validation succeeds.

## First launch

Open **Environment / Database Status**. Capability states come directly from
PhageMine's existing doctor and resource validation logic. A path alone is not
enough for READY: executability, database sidecars, mappings, indexes, and
checksums are validated where configured.

On macOS, Core Analysis uses the bundled PHANOTATE backend. On Windows it stays
NOT READY until a validated backend is configured. Cohort Discovery also needs MMseqs2 for PMF
clustering. Standard Evidence needs ready PHROGs plus MMseqs2. Full Evidence
needs validated Pfam, VOGDB, Swiss-Prot, and PHROGs resources together with
HMMER, DIAMOND, and MMseqs2. `table2asn` is optional.

## Drag-and-drop analysis

Drag one or more supported FASTA files onto the Home upload area. Uploaded names
are reduced to safe portable basenames, FASTA content is validated before run
start, and data is staged in a private temporary directory. Exactly one file is
passed to the single-genome CLI. Two or more files are passed to `phagemine
batch`. Temporary uploaded inputs are retained while the analysis process needs
them and removed after completion.

Filesystem paths are available only under **Advanced input options**. The exact
argument-array command is available under **Advanced / command preview**.

## Evidence database setup

Large databases are never silently bundled or downloaded. No official download
manifest with pinned checksums exists in this repository, so the release
candidate intentionally provides local registration instead of inventing URLs.

Use **Environment / Database Status → Setup / Configure Evidence Databases** to
register a prepared resource and required annotation/metadata sidecar. PhageMine
then runs the existing registry validation and reports READY or specific errors.
Database storage remains outside Applications or Program Files.

## Updating databases

Download updates from the database's official provider, prepare its indexes with
the required upstream tool, and register the new path under a versioned name.
Keep the previous version until representative analyses validate the update.
PhageMine does not overwrite or mutate registered biological databases.

## Troubleshooting

- **Start Analysis is disabled:** read the readiness message and Environment page. Missing dependencies remain missing.
- **FASTA is rejected:** use UTF-8 text with exactly one FASTA record per uploaded file and standard IUPAC DNA symbols.
- **Discovery is unavailable:** install and validate MMseqs2; cohort PMF clustering requires it even with core evidence.
- **Windows full evidence is unavailable:** verify WSL2 and the scientific backend; native HMMER is not claimed.
- **Browser did not open:** launch PhageMine again. The server uses a free loopback port and never exposes the interface to the LAN.
- **Run failed:** open Run Monitor and inspect the diagnostic log and persisted PhageMine checkpoint state. No percentage or result is fabricated.
- **Closing:** in Desktop mode, use **Quit PhageMine** in the sidebar, then close the browser tab.

## CLI for advanced users

The CLI remains fully supported and unchanged. Python users may install the GUI
extra and launch it directly:

```bash
pip install "phagemine[gui]"
phagemine-gui
# or
phagemine gui
```

Scientific commands retain their existing behavior. See the main README for
annotation, discovery, batch, resume, extraction, and database registry syntax.

## Build architecture

`scripts/build_desktop.py` creates native PyInstaller onedir bundles. PyInstaller
is not a cross-compiler, so `.github/workflows/desktop-builds.yml` builds and
smoke-tests Apple Silicon macOS, Intel macOS, and Windows independently. Each
native macOS job wraps `PhageMine.app` in a clearly architecture-specific DMG;
no universal2 claim is made. Windows CI uses Inno Setup to create a Start-menu/desktop installer. Neither
job publishes a GitHub release.

Both native jobs run the loopback-only Streamlit health check and exercise the
frozen CLI against `examples/demo_phage.fasta`. macOS runs a scientific
end-to-end annotation with bundled PHANOTATE 1.6.7 and verifies non-empty native
annotation, evidence, protein/CDS FASTA, GFF3, and manifest outputs. Windows
verifies the same native output pipeline only with the explicitly test-only demo
caller, and separately asserts that doctor does not report Core Analysis READY.
That is packaging evidence, not Windows production gene-calling validation.
