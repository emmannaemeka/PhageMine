# Release checklist

## GitHub source and wheel release

1. Synchronize the stable version in package metadata, citation metadata, Zenodo
   metadata and the local Conda recipe. Add dated release notes and a changelog.
2. Run the complete unit suite and GUI tests on the release commit. Install the
   built wheel with dependencies in clean Linux, Intel macOS and Apple-silicon
   macOS environments; run `pip check` and the installed-package smoke check.
3. Pass Linux external-tool integration and container annotation, output
   validation and bundle verification on that same commit. Preserve check URLs
   and exact dependencies in the release artifacts.
4. Retain the fresh all-six-resource benchmark and sequence-bound naming review.
   Rerun biological searches when functional code or sequence inputs change;
   a metadata-only stable version promotion does not require repeated searches.
5. Build the wheel and source distribution, verify their version metadata and
   SHA256 checksums, and execute the newly built wheel outside source imports.
6. Create the version tag from the checked commit. Upload all assets to a draft
   release, then publish only if the release commit remains current `main` and
   all required checks remain successful. Never reuse a published version.
7. Verify the public release, tag target and downloaded asset checksums.

The publication workflow enforces the current-main CI, external integration and
container gates. Tests use disposable database registries. Release validation
is software validation; it does not establish superior biological accuracy.

## Additional distribution channels

PyPI and Bioconda publication are separate actions. Do not describe them as
available merely because a GitHub release exists. A public Bioconda submission
must replace the local recipe source with the release URL and actual checksum,
then pass channel-specific clean-install and executable checks. Keep external
executables and post-install database setup explicit in every installation route.

## Scientific claims

Use independent representative or held-out review before claiming accuracy
superiority or calibrated confidence. Report reference concordance separately
from functional correctness, and retain unresolved adjudications. Historical
release records remain separate from current evidence.
