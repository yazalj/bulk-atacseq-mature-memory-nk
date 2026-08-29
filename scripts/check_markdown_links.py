#!/usr/bin/env python3
"""Check that repository-local links in Markdown files resolve."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote, urlsplit


ROOT = Path(__file__).resolve().parents[1]
LINK = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
SCHEME = re.compile(r"^[A-Za-z][A-Za-z0-9+.-]*:")


def markdown_files() -> list[Path]:
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
    return sorted(
        Path(line.strip())
        for line in result.stdout.splitlines()
        if line.strip().lower().endswith(".md")
    )


def normalized_target(raw_target: str) -> str | None:
    target = raw_target.strip()
    if target.startswith("<") and ">" in target:
        target = target[1 : target.index(">")]
    elif " \"" in target:
        target = target.split(" \"", maxsplit=1)[0]
    elif " '" in target:
        target = target.split(" '", maxsplit=1)[0]

    if not target or target.startswith(("#", "/", "//")) or SCHEME.match(target):
        return None
    return unquote(urlsplit(target).path)


def main() -> int:
    failures: list[str] = []
    checked = 0
    for relative_source in markdown_files():
        source = ROOT / relative_source
        text = source.read_text(encoding="utf-8")
        for match in LINK.finditer(text):
            target = normalized_target(match.group(1))
            if target is None:
                continue
            checked += 1
            resolved = (source.parent / target).resolve()
            if not resolved.exists():
                line = text.count("\n", 0, match.start()) + 1
                failures.append(f"{relative_source.as_posix()}:{line} -> {target}")

    if failures:
        print("Broken repository-local Markdown links:", file=sys.stderr)
        for failure in failures:
            print(f"- {failure}", file=sys.stderr)
        return 1

    print(f"Markdown link check: PASS ({checked} local links)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
