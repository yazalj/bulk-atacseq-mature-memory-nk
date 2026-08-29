#!/usr/bin/env python3
"""Run the repository's dependency-free synthetic count-matrix smoke test."""

from __future__ import annotations

import argparse
import csv
import json
import math
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SAMPLES = ROOT / "demo" / "input" / "synthetic_samples.tsv"
DEFAULT_COUNTS = ROOT / "demo" / "input" / "synthetic_counts.tsv"
EXPECTED = ROOT / "demo" / "expected"
OUTPUT_NAMES = (
    "retained_counts.tsv",
    "descriptive_summary.tsv",
    "run_summary.json",
)


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Missing header: {path}")
        return reader.fieldnames, list(reader)


def validate_samples(rows: list[dict[str, str]]) -> tuple[list[str], list[str]]:
    required = {"sample_id", "donor", "condition"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("Sample metadata must contain sample_id, donor, and condition.")

    sample_ids = [row["sample_id"] for row in rows]
    if any(not value for value in sample_ids) or len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Sample IDs must be nonempty and unique.")

    conditions = list(dict.fromkeys(row["condition"] for row in rows))
    if len(conditions) != 2:
        raise ValueError("The synthetic demonstration requires exactly two conditions.")

    donors: dict[str, set[str]] = {}
    for row in rows:
        donors.setdefault(row["donor"], set()).add(row["condition"])
    incomplete = [donor for donor, observed in donors.items() if observed != set(conditions)]
    if incomplete:
        raise ValueError(f"Every synthetic donor must have both conditions: {incomplete}")
    return sample_ids, conditions


def parse_counts(
    header: list[str], rows: list[dict[str, str]], sample_ids: list[str]
) -> list[tuple[str, list[int]]]:
    if header != ["peak_id", *sample_ids]:
        raise ValueError("Count-matrix columns must exactly match ordered sample IDs.")

    parsed: list[tuple[str, list[int]]] = []
    seen: set[str] = set()
    for row in rows:
        peak_id = row["peak_id"]
        if not peak_id or peak_id in seen:
            raise ValueError("Peak IDs must be nonempty and unique.")
        seen.add(peak_id)
        try:
            values = [int(row[sample_id]) for sample_id in sample_ids]
        except ValueError as error:
            raise ValueError(f"Counts must be integers for {peak_id}.") from error
        if any(value < 0 for value in values):
            raise ValueError(f"Counts must be nonnegative for {peak_id}.")
        parsed.append((peak_id, values))
    if not parsed:
        raise ValueError("The count matrix contains no peaks.")
    return parsed


def write_tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def run_demo(samples: Path, counts: Path, output_dir: Path) -> dict[str, object]:
    _, sample_rows = read_tsv(samples)
    sample_ids, conditions = validate_samples(sample_rows)
    count_header, count_rows = read_tsv(counts)
    parsed = parse_counts(count_header, count_rows, sample_ids)

    minimum_count = 10
    minimum_samples = 3
    retained = [
        (peak_id, values)
        for peak_id, values in parsed
        if sum(value >= minimum_count for value in values) >= minimum_samples
    ]
    if not retained:
        raise ValueError("The demonstration filter retained no peaks.")

    indexes = {
        condition: [
            index
            for index, row in enumerate(sample_rows)
            if row["condition"] == condition
        ]
        for condition in conditions
    }
    denominator, numerator = conditions
    descriptive_rows: list[list[object]] = []
    for peak_id, values in retained:
        denominator_mean = sum(values[index] for index in indexes[denominator]) / len(
            indexes[denominator]
        )
        numerator_mean = sum(values[index] for index in indexes[numerator]) / len(
            indexes[numerator]
        )
        descriptive_log2fc = math.log2((numerator_mean + 1) / (denominator_mean + 1))
        descriptive_rows.append(
            [
                peak_id,
                f"{denominator_mean:.6f}",
                f"{numerator_mean:.6f}",
                f"{descriptive_log2fc:.6f}",
            ]
        )

    output_dir.mkdir(parents=True, exist_ok=False)
    write_tsv(
        output_dir / "retained_counts.tsv",
        ["peak_id", *sample_ids],
        [[peak_id, *values] for peak_id, values in retained],
    )
    write_tsv(
        output_dir / "descriptive_summary.tsv",
        [
            "peak_id",
            f"mean_{denominator}",
            f"mean_{numerator}",
            f"descriptive_log2fc_{numerator}_vs_{denominator}",
        ],
        descriptive_rows,
    )
    summary: dict[str, object] = {
        "schema_version": 1,
        "samples": len(sample_ids),
        "conditions": conditions,
        "input_peaks": len(parsed),
        "minimum_count": minimum_count,
        "minimum_samples": minimum_samples,
        "retained_peaks": len(retained),
        "interpretation": "technical_smoke_test_only_no_statistical_inference",
    }
    (output_dir / "run_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    return summary


def compare_expected(output_dir: Path) -> None:
    failures = []
    for name in OUTPUT_NAMES:
        observed = (output_dir / name).read_bytes()
        expected = (EXPECTED / name).read_bytes()
        if observed != expected:
            failures.append(name)
    if failures:
        raise ValueError(f"Synthetic outputs differ from expected files: {failures}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=Path, default=DEFAULT_SAMPLES)
    parser.add_argument("--counts", type=Path, default=DEFAULT_COUNTS)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check", action="store_true", help="compare a temporary run with expected outputs")
    mode.add_argument("--output-dir", type=Path, help="write a new demonstration attempt")
    arguments = parser.parse_args()

    if arguments.check:
        with tempfile.TemporaryDirectory(prefix="bulk_atacseq_demo_") as temporary:
            summary = run_demo(arguments.samples, arguments.counts, Path(temporary) / "output")
            compare_expected(Path(temporary) / "output")
    else:
        summary = run_demo(arguments.samples, arguments.counts, arguments.output_dir)

    print(
        f"Synthetic demo: PASS ({summary['retained_peaks']} of "
        f"{summary['input_peaks']} peaks retained)"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
