# PhageMine v1.2 release validation

Validation date: 2026-09-16. A clean temporary clone was created outside the
repository; no benchmark evidence was copied into or changed by the validation
clone.

## Results

- Installed package: **PASS**. Editable installation succeeded in a clean
  project virtual environment with `--no-build-isolation --no-deps`; installed
  version was `1.2.0`.
- Import: **PASS** (`import phagemine`, version `1.2.0`).
- CLI help: **PASS** (`python -m phagemine --help`).
- INPHARED installer dry-run: **PASS**; plan reports the pinned resource and
  approximately 682 MiB compressed download.
- Tests: **PASS**, 272 passed and 10 deselected in the clean clone. The
  in-repository INPHARED/resource test subset passes separately (31 passed,
  including the small-genome regression tests).
- Example annotation execution: intentionally not run because release
  instructions prohibit rerunning annotation tools during this preparation;
  the example input and command are documented for a user environment with
  PHANOTATE and registered databases.

## Environment note

The default isolated pip installation attempt could not fetch build dependency
`setuptools>=68` because this validation environment has no network access.
This is an external network/build-isolation limitation, not a package failure;
the same clean clone installed successfully without build isolation using the
already available build backend.

## Documentation checks

README links to `docs/BENCHMARK.md` and `docs/INPHARED.md`; the benchmark page
links to the distributed compact package and local figure files. Large
third-party databases, raw benchmark evidence, machine-specific logs, and
temporary outputs are excluded by `.gitignore` and remain preserved locally.
