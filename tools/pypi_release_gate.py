#!/usr/bin/env python3
"""Allow PyPI publication only for the published, verified release on main."""
import argparse
import os
import re
from pathlib import Path

from publish_verified_release import api, ready, tag_matches


def published_and_verified(repo, sha, version, get=api):
    if not re.fullmatch(r"[0-9]+\.[0-9]+\.[0-9]+", version):
        return False
    if not ready(repo, sha, get):
        return False
    tag = "v" + version
    release = get(f"repos/{repo}/releases/tags/{tag}", optional=True)
    if not release or release.get("draft", True) or release.get("prerelease", True):
        return False
    if release.get("tag_name") != tag:
        return False
    ref = get(f"repos/{repo}/git/ref/tags/{tag}", optional=True)
    if not ref:
        return False
    try:
        tag_matches(repo, tag, sha, get)
    except ValueError:
        return False
    return True


def main():
    import tomllib
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", required=True)
    parser.add_argument("--sha", required=True)
    args = parser.parse_args()
    version = tomllib.loads(Path("pyproject.toml").read_text())["project"]["version"]
    permitted = published_and_verified(args.repo, args.sha, version)
    if os.environ.get("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"], "a") as handle:
            handle.write("ready=" + str(permitted).lower() + "\n")
    print("READY" if permitted else "SKIPPED_NO_MATCHING_VERIFIED_PUBLISHED_RELEASE")


if __name__ == "__main__":
    main()
