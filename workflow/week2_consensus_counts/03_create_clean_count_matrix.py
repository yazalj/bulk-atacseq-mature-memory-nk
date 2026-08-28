#!/usr/bin/env python3
"""Create an unfiltered clean ATAC-seq count matrix and peak annotation.

Refactored reference version: this is not the script that produced attempt1.
It uses only the Python standard library and performs no normalization,
DESeq2 analysis, or differential-accessibility test.
"""

from __future__ import annotations

import csv
import hashlib
import os
from pathlib import Path
import socket
import sys


PROJECT_ROOT_VALUE = os.environ.get("PROJECT_ROOT")
if not PROJECT_ROOT_VALUE:
    raise SystemExit("Set PROJECT_ROOT to the repository root")
PROJECT_ROOT = Path(PROJECT_ROOT_VALUE)
EXPECTED_HOST = os.environ.get("EXPECTED_HOST", "")
ATTEMPT = os.environ.get("ATTEMPT", "refactored")
EXPECTED_FEATURES = int(os.environ.get("EXPECTED_FEATURES", "112759"))
MANIFEST = Path(
    os.environ.get(
        "MANIFEST",
        str(PROJECT_ROOT / "config/week2_samples.tsv"),
    )
)

ATTEMPT_ROOT = (
    PROJECT_ROOT / "results/week2_consensus_counts" / ATTEMPT
)
SAF = ATTEMPT_ROOT / "annotation/mature_memory_consensus.d0.saf"
RAW = (
    ATTEMPT_ROOT
    / "counts/raw/mature_memory_featurecounts.raw.tsv"
)
SUMMARY = Path(f"{RAW}.summary")
CLEAN_DIR = ATTEMPT_ROOT / "counts/clean"
CLEAN_MATRIX = CLEAN_DIR / "mature_memory_counts.tsv"
ANNOTATION = CLEAN_DIR / "mature_memory_peak_annotation.tsv"
VALIDATION = (
    ATTEMPT_ROOT
    / "counts/validation/clean_count_matrix_validation.tsv"
)
LOG_DIR = (
    PROJECT_ROOT
    / "logs/week2_consensus_counts"
    / ATTEMPT
    / "clean_matrix"
)
CHECKSUMS = LOG_DIR / "03_create_clean_count_matrix_refactored.sha256"

MANIFEST_HEADER = (
    "order",
    "sample_id",
    "donor",
    "cell_type",
    "pairing_status",
    "broadpeak_file",
    "bam_file",
)
SAF_HEADER = ("GeneID", "Chr", "Start", "End", "Strand")
RAW_ANNOTATION_HEADER = ("Geneid", "Chr", "Start", "End", "Strand", "Length")


class ValidationError(RuntimeError):
    """Raised when an input or output does not match the locked workflow."""


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require_source(path: Path, label: str) -> None:
    if path.is_symlink() or not path.is_file() or path.stat().st_size == 0:
        raise ValidationError(f"{label} must be a nonempty regular file: {path}")


def read_manifest() -> tuple[list[str], list[str]]:
    with MANIFEST.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if tuple(reader.fieldnames or ()) != MANIFEST_HEADER:
            raise ValidationError("sample-manifest header is not exact")
        rows = list(reader)

    if len(rows) != 10:
        raise ValidationError("sample manifest must contain exactly ten rows")
    if [row["order"] for row in rows] != [str(i) for i in range(1, 11)]:
        raise ValidationError("sample-manifest order must be 1 through 10")

    sample_ids = [row["sample_id"] for row in rows]
    bam_names = [row["bam_file"] for row in rows]
    if len(set(sample_ids)) != 10 or len(set(bam_names)) != 10:
        raise ValidationError("sample IDs and BAM basenames must be unique")
    if any("/" in name for name in bam_names):
        raise ValidationError("BAM entries in the manifest must be basenames")

    cell_types = [row["cell_type"] for row in rows]
    if cell_types.count("Mature_NK") != 4 or cell_types.count("Memory_NK") != 6:
        raise ValidationError("manifest must contain four Mature and six Memory samples")
    return sample_ids, bam_names


def read_assigned_counts(expected_bams: list[str]) -> list[int]:
    with SUMMARY.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None or header[0] != "Status":
            raise ValidationError("invalid featureCounts summary header")
        if [Path(value).name for value in header[1:]] != expected_bams:
            raise ValidationError("summary BAM columns do not match the locked order")

        statuses: set[str] = set()
        assigned: list[int] | None = None
        categories = 0
        for row in reader:
            categories += 1
            if len(row) != 11 or not row[0] or row[0] in statuses:
                raise ValidationError("invalid or duplicate summary row")
            statuses.add(row[0])
            if any(not value.isascii() or not value.isdigit() for value in row[1:]):
                raise ValidationError("summary contains a non-integer count")
            if row[0] == "Assigned":
                assigned = [int(value) for value in row[1:]]

    if categories != 14 or assigned is None:
        raise ValidationError("summary must contain 14 categories and one Assigned row")
    return assigned


def raw_header(reader: csv.reader) -> list[str]:
    comments = 0
    for row in reader:
        if row and row[0].startswith("#"):
            comments += 1
            continue
        if comments < 1 or not row:
            raise ValidationError("raw featureCounts header is missing")
        return row
    raise ValidationError("raw featureCounts file has no header")


def write_outputs(
    sample_ids: list[str],
    expected_bams: list[str],
    assigned: list[int],
) -> tuple[int, list[int]]:
    with (
        RAW.open("r", encoding="utf-8", newline="") as raw_handle,
        SAF.open("r", encoding="utf-8", newline="") as saf_handle,
        CLEAN_MATRIX.open("x", encoding="utf-8", newline="") as clean_handle,
        ANNOTATION.open("x", encoding="utf-8", newline="") as annotation_handle,
    ):
        raw_reader = csv.reader(raw_handle, delimiter="\t")
        saf_reader = csv.reader(saf_handle, delimiter="\t")
        raw_columns = raw_header(raw_reader)
        if tuple(raw_columns[:6]) != RAW_ANNOTATION_HEADER:
            raise ValidationError("raw annotation columns are not exact")
        if [Path(value).name for value in raw_columns[6:]] != expected_bams:
            raise ValidationError("raw BAM columns do not match the locked order")
        if tuple(next(saf_reader, ())) != SAF_HEADER:
            raise ValidationError("SAF header is not exact")

        clean_writer = csv.writer(
            clean_handle,
            delimiter="\t",
            lineterminator="\n",
        )
        annotation_writer = csv.writer(
            annotation_handle,
            delimiter="\t",
            lineterminator="\n",
        )
        clean_writer.writerow(("peak_id", *sample_ids))
        annotation_writer.writerow(
            ("peak_id", "Chr", "Start", "End", "Strand", "Length")
        )

        matrix_sums = [0] * 10
        seen_ids: set[str] = set()
        rows = 0

        for raw_row in raw_reader:
            rows += 1
            saf_row = next(saf_reader, None)
            if saf_row is None or len(raw_row) != 16 or len(saf_row) != 5:
                raise ValidationError(f"raw/SAF structure mismatch at row {rows}")

            peak_id, chromosome, start, end, strand, length = raw_row[:6]
            expected_id = f"peak_{rows:06d}"
            if peak_id != expected_id or peak_id in seen_ids:
                raise ValidationError(f"invalid or duplicate peak ID at row {rows}")
            seen_ids.add(peak_id)
            if saf_row != [peak_id, chromosome, start, end, strand]:
                raise ValidationError(f"raw annotation differs from SAF at row {rows}")
            if (
                not start.isdigit()
                or not end.isdigit()
                or not length.isdigit()
                or int(start) < 1
                or int(end) < int(start)
                or int(length) != int(end) - int(start) + 1
                or strand != "."
            ):
                raise ValidationError(f"invalid annotation at row {rows}")

            count_text = raw_row[6:]
            if any(not value.isascii() or not value.isdigit() for value in count_text):
                raise ValidationError(f"non-integer count at row {rows}")
            counts = [int(value) for value in count_text]
            for index, value in enumerate(counts):
                matrix_sums[index] += value

            clean_writer.writerow((peak_id, *count_text))
            annotation_writer.writerow(
                (peak_id, chromosome, start, end, strand, length)
            )

        if next(saf_reader, None) is not None:
            raise ValidationError("SAF contains extra rows")
        if rows != EXPECTED_FEATURES:
            raise ValidationError(
                f"observed {rows} features; expected {EXPECTED_FEATURES}"
            )
        if matrix_sums != assigned:
            raise ValidationError("clean matrix sums do not match summary Assigned")
    return rows, matrix_sums


def main() -> int:
    if not ATTEMPT or any(character not in "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789._-" for character in ATTEMPT):
        raise ValidationError(f"unsafe ATTEMPT name: {ATTEMPT!r}")
    if EXPECTED_HOST and socket.gethostname().split(".")[0] != EXPECTED_HOST:
        raise ValidationError(
            f"expected host {EXPECTED_HOST}; observed {socket.gethostname()}"
        )
    if Path.cwd().resolve() != PROJECT_ROOT.resolve():
        raise ValidationError(f"run from exactly {PROJECT_ROOT}")

    for path, label in (
        (MANIFEST, "sample manifest"),
        (SAF, "consensus SAF"),
        (RAW, "raw featureCounts table"),
        (SUMMARY, "featureCounts summary"),
    ):
        require_source(path, label)

    if CLEAN_DIR.exists() or CLEAN_DIR.is_symlink():
        raise ValidationError(f"refusing to overwrite clean directory: {CLEAN_DIR}")
    if LOG_DIR.exists() or LOG_DIR.is_symlink():
        raise ValidationError(f"refusing to overwrite log directory: {LOG_DIR}")
    if VALIDATION.exists() or VALIDATION.is_symlink():
        raise ValidationError(f"refusing to overwrite validation file: {VALIDATION}")

    source_paths = (MANIFEST, SAF, RAW, SUMMARY)
    source_hashes = {path: sha256(path) for path in source_paths}
    sample_ids, bam_names = read_manifest()
    assigned = read_assigned_counts(bam_names)

    CLEAN_DIR.mkdir()
    LOG_DIR.mkdir(parents=True)
    rows, matrix_sums = write_outputs(sample_ids, bam_names, assigned)

    if any(sha256(path) != source_hashes[path] for path in source_paths):
        raise ValidationError("one or more source files changed during conversion")

    with VALIDATION.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("check", "value", "status"))
        writer.writerow(("feature_rows", rows, "PASS"))
        writer.writerow(("sample_count", len(sample_ids), "PASS"))
        writer.writerow(("peak_and_sample_order_preserved", "YES", "PASS"))
        writer.writerow(("counts_match_raw_and_Assigned", "YES", "PASS"))
        writer.writerow(("missing_count_values", 0, "PASS"))
        writer.writerow(("filtering_applied", "NO", "PASS"))
        writer.writerow(("normalization_applied", "NO", "PASS"))
        for sample_id, value in zip(sample_ids, matrix_sums):
            writer.writerow((f"matrix_sum_{sample_id}", value, "PASS"))

    checksum_targets = (CLEAN_MATRIX, ANNOTATION, VALIDATION)
    with CHECKSUMS.open("x", encoding="utf-8") as handle:
        for path in checksum_targets:
            handle.write(f"{sha256(path)}  {path}\n")

    print(f"clean_matrix={CLEAN_MATRIX}")
    print(f"annotation={ANNOTATION}")
    print(f"features={rows}")
    print("filtering_applied=NO")
    print("normalization_applied=NO")
    print("DESeq2_executed=NO")
    print("status=PASS")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, ValidationError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
