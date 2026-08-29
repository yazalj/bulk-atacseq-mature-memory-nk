#!/usr/bin/env python3
"""Write or verify the repository SHA-256 manifest."""

from __future__ import annotations

import argparse
import hashlib
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "MANIFEST.sha256"


def repository_files() -> list[Path]:
    command = ["git", "-c", f"safe.directory={ROOT.as_posix()}", "ls-files", "--cached"]
    result = subprocess.run(
        command,
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    paths = {
        Path(line.strip())
        for line in result.stdout.splitlines()
        if line.strip() and line.strip() != MANIFEST.name
    }
    return sorted(paths, key=lambda path: path.as_posix())


def tracked_sha256(path: Path) -> str:
    result = subprocess.run(
        [
            "git",
            "-c",
            f"safe.directory={ROOT.as_posix()}",
            "show",
            f":{path.as_posix()}",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
    )
    return hashlib.sha256(result.stdout).hexdigest()


def write_manifest() -> int:
    paths = repository_files()
    lines = [f"{tracked_sha256(path)}  {path.as_posix()}" for path in paths]
    MANIFEST.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"Wrote {len(lines)} checksums to {MANIFEST.name}.")
    return 0


def read_manifest() -> dict[Path, str]:
    entries: dict[Path, str] = {}
    for number, raw_line in enumerate(
        MANIFEST.read_text(encoding="utf-8").splitlines(), start=1
    ):
        if not raw_line.strip():
            continue
        try:
            digest, name = raw_line.split("  ", maxsplit=1)
        except ValueError as error:
            raise ValueError(f"Malformed manifest line {number}: {raw_line}") from error
        if len(digest) != 64 or any(char not in "0123456789abcdef" for char in digest):
            raise ValueError(f"Invalid SHA-256 on manifest line {number}")
        path = Path(name)
        if path in entries:
            raise ValueError(f"Duplicate manifest path: {name}")
        entries[path] = digest
    return entries


def check_manifest() -> int:
    expected_paths = set(repository_files())
    entries = read_manifest()
    listed_paths = set(entries)
    failures: list[str] = []

    for path in sorted(expected_paths - listed_paths, key=lambda item: item.as_posix()):
        failures.append(f"not listed: {path.as_posix()}")
    for path in sorted(listed_paths - expected_paths, key=lambda item: item.as_posix()):
        failures.append(f"not tracked: {path.as_posix()}")
    for path in sorted(expected_paths & listed_paths, key=lambda item: item.as_posix()):
        full_path = ROOT / path
        if not full_path.is_file():
            failures.append(f"missing: {path.as_posix()}")
            continue
        actual = tracked_sha256(path)
        if actual != entries[path]:
            failures.append(f"checksum mismatch: {path.as_posix()}")

    if failures:
        print("Manifest verification failed:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1

    print(f"Manifest verification: PASS ({len(entries)} files)")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true", help="verify the manifest")
    mode.add_argument("--write", action="store_true", help="regenerate the manifest")
    arguments = parser.parse_args()
    return write_manifest() if arguments.write else check_manifest()


if __name__ == "__main__":
    raise SystemExit(main())
