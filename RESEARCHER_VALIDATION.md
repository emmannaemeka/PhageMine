# Researcher usability and software validation

Date: 2026-10-03. Branch: `improve-researcher-readiness`.
Base commit: `b407a0c32eea060fab5d5d7ce2c69769a515771c`.
These are development-branch checks, not a replacement for the historical
v1.2.0 release record or evidence of biological superiority.

## Verified locally

| Check | Result | Scope |
|---|---|---|
| Fresh environment installation with test and GUI dependencies | PASS | Normal isolated pip build and dependency resolution on Linux, Python 3.12.14 |
| Baseline unit suite | PASS | 322 passed, 2 skipped, 12 integration tests deselected |
| Updated unit suite | PASS | 336 passed, 2 skipped, 12 integration tests deselected |
| Integration-selected suite | PASS with incomplete coverage | 11 passed, 1 skipped because its external resource was unavailable |
| Source and wheel build | PASS | Built through the standard isolated build backend |
| Built wheel installed in a second fresh environment | PASS | All declared runtime and GUI dependencies resolved; pip check found no broken requirements |
| Installed synthetic workflow outside checkout | PASS | Five records; 15 nonempty outputs and consistent annotation/evidence/classification IDs |
| Installed real-PHANOTATE Core workflow outside checkout | PASS | Three records on the synthetic bundled input; same output-integrity checks |
| GUI service and application tests | PASS | Native layouts, status reporting, file validation, exports, page rendering and successful-run UI |
| Actual browser | PASS for inspected pages | Home, native annotation details, saved status, figures and export pages; screenshots retained |

The smoke-check script uses a temporary working directory, removes PYTHONPATH,
and supplies a disposable resource registry. No mock evidence is enabled.
Synthetic-predictor and real-caller checks are labelled separately and both
explicitly state that scientific accuracy was not validated.

## Installation failure found and corrected

PHANOTATE 1.6.7 failed to import `pkg_resources` in a fresh Python 3.12
environment. Installing `setuptools<81` fixed the import and the Core run
completed. That compatibility requirement is now explicit in environment.yml
and installation troubleshooting. The pip installation of PHANOTATE's
fastpath dependency also required a compiler; local validation used GCC.
The documented Conda route supplies the external caller as a packaged tool.

## Not established by these checks

- Clean Conda installation on Linux and both macOS architectures was not run
  locally. Remote CI has been strengthened but has not been executed here.
- The six large evidence databases were not installed. Full-evidence annotation,
  real discovery/Both cohorts, resource-backed resume and scale checks remain
  release gates. Passing tests with fixtures does not certify those deployments.
- Independent expert adjudication, leakage-controlled held-out evaluation,
  calibrated confidence and universal superiority remain outstanding.
- The old seven-genome results and the 476/689 unresolved adjudications were
  preserved; no validation result was fabricated or relabelled.
- The GUI still waits for the analysis subprocess and constructs ZIPs in memory
  on request. Background job control and large-result export remain limitations.
- The published v1.2.0 artifact and Zenodo record have not been changed. A new
  release needs its own version and completed release gates.

Generated installation logs, distributions and run outputs are under the
ignored `.validation/` directory. Development screenshots are in
`docs/images/phagemine-gui/`. See [the readiness assessment](docs/researcher-readiness.md)
and [release checklist](docs/release-checklist.md) for outstanding work.
