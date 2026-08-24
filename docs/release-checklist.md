# Release checklist

A PhageMine release is complete only after the published artifact, not merely
the source checkout, passes these checks.

1. Update the version in `pyproject.toml`, `src/phagemine/__init__.py`,
   `CITATION.cff`, and `CHANGELOG.md`; run the metadata consistency test.
2. Run the complete unit suite and the external-tool integration suite.
3. Create the release tag from the exact reviewed commit. Do not add documented
   functionality to `main` under an already published version number.
4. Update the Bioconda recipe to the new tag and checksum. Its run dependencies
   must include `phanotate`, `hmmer`, `mmseqs2`, `diamond`, and `mash`.
5. In clean Linux, Intel macOS, and Apple-silicon macOS environments, install
   the Bioconda artifact using the documented command.
6. Run each binary directly and require a zero exit status: PHANOTATE, HMMER,
   MMseqs2, DIAMOND, and Mash.
7. Run `phagemine doctor`; all production executables and intended capabilities
   must be `READY`, and the command must exit zero.
8. Install all six databases in a clean user-data directory and run Doctor
   again. Repeat with a registry retained from the immediately previous
   PhageMine version.
9. Run the bundled demo genome through Core and full-evidence annotation, and
   run a small cohort through discovery/Both mode. Verify expected reports,
   tables, figures, provenance, and GenBank pre-submission outputs.
10. Publish the release only after the clean-install evidence and exact package
    versions are recorded in the release notes.
