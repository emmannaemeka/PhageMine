#!/usr/bin/env python3
"""Fail-closed architecture and deployment-target audit for a macOS app."""
from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path


def command(*args: str) -> str:
    return subprocess.run(args, check=True, capture_output=True, text=True).stdout


def is_macho(path: Path) -> bool:
    return "Mach-O" in command("file", "-b", str(path))


def architectures(path: Path) -> set[str]:
    return set(command("lipo", "-archs", str(path)).strip().split())


def parse_minimum_versions(load_commands: str) -> list[tuple[int, ...]]:
    """Read only deployment-target load commands, not dylib version fields."""
    values: list[str] = []
    command_name: str | None = None
    for line in load_commands.splitlines():
        match = re.match(r"^\s*cmd\s+(\S+)\s*$", line)
        if match:
            command_name = match.group(1)
            continue
        if command_name == "LC_BUILD_VERSION":
            match = re.match(r"^\s*minos\s+(\d+(?:\.\d+){1,2})\s*$", line)
        elif command_name == "LC_VERSION_MIN_MACOSX":
            match = re.match(r"^\s*version\s+(\d+(?:\.\d+){1,2})\s*$", line)
        else:
            match = None
        if match:
            values.append(match.group(1))
    return [tuple(int(part) for part in value.split(".")) for value in values]


def minimum_versions(path: Path) -> list[tuple[int, ...]]:
    return parse_minimum_versions(command("otool", "-l", str(path)))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("bundle", type=Path)
    parser.add_argument("--architecture", choices=("arm64", "x86_64"), required=True)
    parser.add_argument("--maximum-minimum-macos")
    args = parser.parse_args()

    required = [
        args.bundle / "Contents" / "MacOS" / "PhageMine",
        args.bundle / "Contents" / "Resources" / "phagemine_backend" / "PhageMine-PHANOTATE",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise SystemExit(f"Required packaged executables are missing: {missing}")

    macho_files = [path for path in args.bundle.rglob("*") if path.is_file() and is_macho(path)]
    if not macho_files:
        raise SystemExit("No Mach-O files found in application bundle")

    maximum = tuple(int(part) for part in args.maximum_minimum_macos.split(".")) if args.maximum_minimum_macos else None
    errors: list[str] = []
    for path in macho_files:
        found_architectures = architectures(path)
        if args.architecture not in found_architectures:
            errors.append(f"{path}: architectures {sorted(found_architectures)!r} omit {args.architecture}")
        if maximum is not None:
            versions = minimum_versions(path)
            if not versions:
                errors.append(f"{path}: no macOS minimum deployment version load command")
            for version in versions:
                padded_version = version + (0,) * (len(maximum) - len(version))
                padded_maximum = maximum + (0,) * (len(version) - len(maximum))
                if padded_version > padded_maximum:
                    errors.append(
                        f"{path}: minimum macOS {'.'.join(map(str, version))} exceeds "
                        f"{args.maximum_minimum_macos}"
                    )
    if errors:
        raise SystemExit("macOS bundle audit failed:\n" + "\n".join(errors))

    print(
        f"Validated {len(macho_files)} Mach-O files for {args.architecture}"
        + (f" with minimum macOS <= {args.maximum_minimum_macos}" if maximum else "")
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
