#!/usr/bin/env python3
"""Stage tracked public files for the MkDocs documentation build."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
STAGE = (ROOT / ".site_source").resolve()
INCLUDED_ROOT_FILES = {
    "AGENTS.md",
    "CHANGELOG.md",
    "CITATION.cff",
    "LICENSE",
    "README.md",
    "THIRD_PARTY_NOTICES.md",
}
INCLUDED_PREFIXES = (
    ".agents/",
    "agent/",
    "config/",
    "demo/",
    "docs/",
    "environment/",
    "results/",
    "scripts/",
    "workflow/",
)


def repository_files() -> list[Path]:
    result = subprocess.run(
        [
            "git",
            "-c",
            f"safe.directory={ROOT.as_posix()}",
            "ls-files",
            "--cached",
            "--others",
            "--exclude-standard",
        ],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    return [Path(line.strip()) for line in result.stdout.splitlines() if line.strip()]


def included(path: Path) -> bool:
    posix = path.as_posix()
    return posix in INCLUDED_ROOT_FILES or posix.startswith(INCLUDED_PREFIXES)


def main() -> int:
    if STAGE.parent != ROOT or STAGE.name != ".site_source":
        raise RuntimeError(f"Refusing to replace unexpected staging path: {STAGE}")
    if STAGE.exists():
        shutil.rmtree(STAGE)
    STAGE.mkdir()

    copied = 0
    for relative_path in repository_files():
        if not included(relative_path):
            continue
        source = ROOT / relative_path
        if not source.is_file():
            continue
        destination_path = Path("index.md") if relative_path.as_posix() == "README.md" else relative_path
        destination = STAGE / destination_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)
        copied += 1

    print(f"Prepared {copied} files in {STAGE.name}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
