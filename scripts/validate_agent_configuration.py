#!/usr/bin/env python3
"""Validate the machine-readable inputs for an AI-assisted ATAC-seq run.

The validator uses only the Python standard library. It is deliberately
read-only: it reports structural errors and review warnings but never changes
configuration, data, analysis state, or outputs.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path
import re
import sys
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CONTRACT = REPO_ROOT / "agent" / "workflow_contract.yaml"

STAGE_IDS = [
    "define_study",
    "validate_inputs",
    "raw_read_qc",
    "trim_reads",
    "align_reads",
    "process_bam",
    "call_peaks",
    "build_consensus",
    "count_fragments",
    "differential_accessibility",
    "downstream_analysis",
    "report_analysis",
]

DECISION_IDS = [
    "study_scope",
    "reference_resources",
    "read_processing",
    "alignment_filtering",
    "peak_calling",
    "consensus_strategy",
    "statistical_design",
    "downstream_scope",
]

VALID_MODES = {
    "case_study_audit",
    "new_study_dry_run",
    "new_study_execution",
}

VALID_STAGE_STATUSES = {
    "not_started",
    "awaiting_human",
    "approved",
    "running",
    "validation_failed",
    "complete",
    "blocked",
    "not_applicable",
}

DEPENDENCIES = {
    "define_study": [],
    "validate_inputs": ["define_study"],
    "raw_read_qc": ["validate_inputs"],
    "trim_reads": ["raw_read_qc"],
    "align_reads": ["raw_read_qc"],
    "process_bam": ["align_reads"],
    "call_peaks": ["process_bam"],
    "build_consensus": ["call_peaks"],
    "count_fragments": ["build_consensus", "process_bam"],
    "differential_accessibility": ["count_fragments"],
    "downstream_analysis": ["differential_accessibility"],
    "report_analysis": ["differential_accessibility"],
}

STAGE_DECISIONS = {
    "define_study": ["study_scope"],
    "validate_inputs": ["reference_resources"],
    "raw_read_qc": ["read_processing"],
    "trim_reads": ["read_processing"],
    "align_reads": ["reference_resources", "alignment_filtering"],
    "process_bam": ["alignment_filtering"],
    "call_peaks": ["peak_calling"],
    "build_consensus": ["consensus_strategy"],
    "count_fragments": ["consensus_strategy"],
    "differential_accessibility": ["statistical_design"],
    "downstream_analysis": ["downstream_scope"],
    "report_analysis": [],
}

SAMPLE_COLUMNS = [
    "order",
    "sample_id",
    "donor",
    "condition",
    "batch",
    "treatment",
    "sequencing_layout",
    "read1_fastq",
    "read2_fastq",
]

PLACEHOLDER = re.compile(
    r"(^|[^a-z])(review|replace)[_-]?(me|from|with)?|"
    r"/path/to/|record_sha256",
    re.IGNORECASE,
)


class Report:
    """Collect deterministic validation messages."""

    def __init__(self) -> None:
        self.errors: list[str] = []
        self.warnings: list[str] = []
        self.passes: list[str] = []

    def error(self, message: str) -> None:
        self.errors.append(message)

    def warn(self, message: str) -> None:
        self.warnings.append(message)

    def passed(self, message: str) -> None:
        self.passes.append(message)


def load_json(path: Path, label: str, report: Report) -> Any | None:
    if not path.is_file():
        report.error(f"{label} does not exist as a regular file: {path}")
        return None
    try:
        with path.open("r", encoding="utf-8") as handle:
            return json.load(handle)
    except (OSError, UnicodeError, json.JSONDecodeError) as exc:
        report.error(f"{label} is not valid UTF-8 JSON: {path}: {exc}")
        return None


def find_placeholders(value: Any, location: str = "root") -> list[str]:
    found: list[str] = []
    if isinstance(value, dict):
        for key, child in value.items():
            found.extend(find_placeholders(child, f"{location}.{key}"))
    elif isinstance(value, list):
        for index, child in enumerate(value):
            found.extend(find_placeholders(child, f"{location}[{index}]"))
    elif isinstance(value, str) and PLACEHOLDER.search(value):
        found.append(location)
    return found


def nested(
    data: Any,
    keys: tuple[str, ...],
    report: Report,
) -> Any | None:
    cursor = data
    for key in keys:
        if not isinstance(cursor, dict) or key not in cursor:
            report.error(f"parameters missing required field: {'.'.join(keys)}")
            return None
        cursor = cursor[key]
    return cursor


def validate_contract(path: Path, report: Report) -> None:
    if not path.is_file():
        report.error(f"workflow contract is missing: {path}")
        return
    text = path.read_text(encoding="utf-8")
    if "\nstages:\n" not in text:
        report.error("workflow contract has no top-level stages section")
        return
    stages_text = text.split("\nstages:\n", maxsplit=1)[1]
    observed = re.findall(
        r"^  - id: ([a-z0-9_]+)$",
        stages_text,
        re.MULTILINE,
    )
    if observed != STAGE_IDS:
        report.error(
            "workflow contract stage order differs from the validator: "
            f"expected {STAGE_IDS}; observed {observed}"
        )
        return
    report.passed("workflow contract contains the expected ordered stage IDs")


def validate_samples(
    path: Path,
    report: Report,
    check_files: bool,
) -> tuple[list[dict[str, str]], set[str]]:
    if not path.is_file():
        report.error(f"sample registry does not exist: {path}")
        return [], set()

    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            header = reader.fieldnames or []
            missing = [column for column in SAMPLE_COLUMNS if column not in header]
            if missing:
                report.error(
                    "sample registry is missing required columns: "
                    + ", ".join(missing)
                )
                return [], set()
            rows = list(reader)
    except (OSError, UnicodeError, csv.Error) as exc:
        report.error(f"sample registry could not be read: {path}: {exc}")
        return [], set()

    if len(rows) < 2:
        report.error("sample registry must contain at least two libraries")
        return rows, set()

    seen_samples: set[str] = set()
    seen_paths: set[str] = set()
    conditions: list[str] = []
    orders: list[int] = []

    for line_number, row in enumerate(rows, start=2):
        for column in SAMPLE_COLUMNS:
            if not (row.get(column) or "").strip():
                report.error(f"samples line {line_number}: empty {column}")

        try:
            order = int((row.get("order") or "").strip())
            orders.append(order)
        except ValueError:
            report.error(f"samples line {line_number}: order is not an integer")

        sample_id = (row.get("sample_id") or "").strip()
        if sample_id in seen_samples:
            report.error(f"samples line {line_number}: duplicate sample_id {sample_id}")
        seen_samples.add(sample_id)

        if (row.get("sequencing_layout") or "").strip() != "paired-end":
            report.error(
                f"samples line {line_number}: sequencing_layout must be paired-end"
            )

        read1 = (row.get("read1_fastq") or "").strip()
        read2 = (row.get("read2_fastq") or "").strip()
        if read1 == read2 and read1:
            report.error(
                f"samples line {line_number}: read1_fastq and read2_fastq are identical"
            )
        for mate, recorded_path in (("R1", read1), ("R2", read2)):
            if recorded_path in seen_paths and recorded_path:
                report.error(
                    f"samples line {line_number}: {mate} path is reused: {recorded_path}"
                )
            seen_paths.add(recorded_path)
            if PLACEHOLDER.search(recorded_path):
                report.error(
                    f"samples line {line_number}: unresolved placeholder in {mate} path"
                )
            if check_files and recorded_path:
                candidate = Path(recorded_path)
                if not candidate.is_absolute():
                    candidate = (Path.cwd() / candidate).resolve()
                if not candidate.is_file() or candidate.stat().st_size == 0:
                    report.error(
                        f"samples line {line_number}: {mate} is missing or empty: "
                        f"{candidate}"
                    )

        conditions.append((row.get("condition") or "").strip())

    expected_orders = list(range(1, len(rows) + 1))
    if orders != expected_orders:
        report.error(
            "sample order must be consecutive and match row order: "
            f"expected {expected_orders}; observed {orders}"
        )

    condition_counts = Counter(conditions)
    condition_set = {condition for condition in conditions if condition}
    if len(condition_set) < 2:
        report.error("at least two conditions are required for a contrast")
    for condition, count in sorted(condition_counts.items()):
        if condition and count < 2:
            report.warn(
                f"condition {condition!r} has fewer than two libraries; "
                "biological replication requires human review"
            )

    if not report.errors:
        report.passed(
            f"sample registry contains {len(rows)} unique paired-end libraries"
        )
    return rows, condition_set


def validate_parameters(
    path: Path,
    rows: list[dict[str, str]],
    conditions: set[str],
    report: Report,
) -> dict[str, Any] | None:
    data = load_json(path, "parameter document", report)
    if not isinstance(data, dict):
        if data is not None:
            report.error("parameter document root must be a JSON object")
        return None

    placeholders = find_placeholders(data)
    for location in placeholders:
        report.error(f"unresolved parameter placeholder at {location}")

    if data.get("schema_version") != 1:
        report.error("parameters.schema_version must equal 1")
    analysis_id = data.get("analysis_id")
    if not isinstance(analysis_id, str) or not analysis_id.strip():
        report.error("parameters.analysis_id must be a nonempty string")

    assay = nested(data, ("project", "assay"), report)
    layout = nested(data, ("project", "sequencing_layout"), report)
    organism = nested(data, ("project", "organism"), report)
    genome_build = nested(data, ("project", "genome_build"), report)
    if assay != "bulk_ATAC_seq":
        report.error("project.assay must be bulk_ATAC_seq")
    if layout != "paired_end":
        report.error("project.sequencing_layout must be paired_end")
    if not isinstance(organism, str) or not organism.strip():
        report.error("project.organism must be a nonempty string")
    if not isinstance(genome_build, str) or not genome_build.strip():
        report.error("project.genome_build must be a nonempty string")

    threads = nested(data, ("compute", "maximum_total_threads"), report)
    concurrent = nested(data, ("compute", "concurrent_intensive_stages"), report)
    authorized_by = nested(data, ("compute", "thread_limit_authorized_by"), report)
    if not isinstance(threads, int) or isinstance(threads, bool) or threads < 1:
        report.error("compute.maximum_total_threads must be a positive integer")
    elif threads > 2:
        if not isinstance(authorized_by, str) or not authorized_by.strip():
            report.error("a thread limit above 2 requires thread_limit_authorized_by")
        else:
            report.warn(
                f"maximum_total_threads is {threads}; confirm the recorded authorization"
            )
    if concurrent != 1:
        report.error("compute.concurrent_intensive_stages must equal 1")

    for key in (
        "fasta",
        "aligner_index_prefix",
        "annotation_gtf",
        "blacklist_bed",
        "chromosome_style",
    ):
        nested(data, ("references", key), report)
    chromosome_style = nested(data, ("references", "chromosome_style"), report)
    if chromosome_style not in {"chr_prefixed", "unprefixed", "other"}:
        report.error(
            "references.chromosome_style must be chr_prefixed, unprefixed, or other"
        )

    trimming = nested(data, ("preprocessing", "perform_trimming"), report)
    adapter1 = nested(data, ("preprocessing", "adapter_read1"), report)
    adapter2 = nested(data, ("preprocessing", "adapter_read2"), report)
    mapq = nested(data, ("preprocessing", "minimum_mapping_quality"), report)
    proper_pair = nested(data, ("preprocessing", "require_proper_pair"), report)
    duplicate_policy = nested(data, ("preprocessing", "duplicate_policy"), report)
    blacklist_policy = nested(data, ("preprocessing", "blacklist_removal"), report)
    blacklist = nested(data, ("references", "blacklist_bed"), report)

    if not isinstance(trimming, bool):
        report.error("preprocessing.perform_trimming must be boolean")
    elif trimming:
        if not isinstance(adapter1, str) or not adapter1.strip():
            report.error("adapter_read1 is required when trimming is enabled")
        if not isinstance(adapter2, str) or not adapter2.strip():
            report.error("adapter_read2 is required when trimming is enabled")
    if not isinstance(mapq, int) or isinstance(mapq, bool) or not 0 <= mapq <= 60:
        report.error("minimum_mapping_quality must be an integer from 0 to 60")
    if proper_pair is not True:
        report.error("require_proper_pair must be true for this paired-end contract")
    if duplicate_policy not in {"mark", "remove", "retain"}:
        report.error("duplicate_policy must be mark, remove, or retain")
    if blacklist_policy not in {"pair_preserving", "not_used"}:
        report.error("blacklist_removal must be pair_preserving or not_used")
    if blacklist_policy == "pair_preserving" and not isinstance(blacklist, str):
        report.error("pair-preserving blacklist removal requires blacklist_bed")
    if blacklist_policy == "not_used" and blacklist is not None:
        report.warn("blacklist_bed is present although blacklist_removal is not_used")

    caller = nested(data, ("peak_calling", "caller"), report)
    peak_format = nested(data, ("peak_calling", "input_format"), report)
    peak_type = nested(data, ("peak_calling", "peak_type"), report)
    peak_qvalue = nested(data, ("peak_calling", "qvalue"), report)
    if caller != "MACS3":
        report.error("peak_calling.caller must be MACS3 for this contract")
    if peak_format != "BAMPE":
        report.error("peak_calling.input_format must be BAMPE")
    if peak_type not in {"narrow", "broad"}:
        report.error("peak_calling.peak_type must be narrow or broad")
    if not isinstance(peak_qvalue, (int, float)) or isinstance(peak_qvalue, bool):
        report.error("peak_calling.qvalue must be numeric")
    elif not 0 < float(peak_qvalue) <= 1:
        report.error("peak_calling.qvalue must be greater than 0 and at most 1")

    count_unit = nested(data, ("counting", "count_unit"), report)
    strandedness = nested(data, ("counting", "strandedness"), report)
    consensus = nested(data, ("counting", "consensus_strategy"), report)
    if count_unit != "paired_fragments":
        report.error("counting.count_unit must be paired_fragments")
    if strandedness != "unstranded":
        report.error("counting.strandedness must be unstranded")
    if not isinstance(consensus, str) or not consensus.strip():
        report.error("counting.consensus_strategy must be a nonempty string")

    minimum_count = nested(
        data, ("differential_accessibility", "minimum_count"), report
    )
    minimum_samples = nested(
        data, ("differential_accessibility", "minimum_samples"), report
    )
    design = nested(data, ("differential_accessibility", "design"), report)
    reference = nested(
        data, ("differential_accessibility", "reference_level"), report
    )
    numerator = nested(
        data, ("differential_accessibility", "contrast_numerator"), report
    )
    denominator = nested(
        data, ("differential_accessibility", "contrast_denominator"), report
    )
    alpha = nested(
        data,
        ("differential_accessibility", "adjusted_pvalue_threshold"),
        report,
    )
    lfc_threshold = nested(
        data,
        ("differential_accessibility", "absolute_log2_fold_change_threshold"),
        report,
    )
    shrinkage = nested(
        data, ("differential_accessibility", "lfc_shrinkage"), report
    )

    if not isinstance(minimum_count, int) or isinstance(minimum_count, bool):
        report.error("minimum_count must be a nonnegative integer")
    elif minimum_count < 0:
        report.error("minimum_count must be nonnegative")
    if not isinstance(minimum_samples, int) or isinstance(minimum_samples, bool):
        report.error("minimum_samples must be a positive integer")
    elif not 1 <= minimum_samples <= len(rows):
        report.error(
            f"minimum_samples must be between 1 and the sample count ({len(rows)})"
        )
    if not isinstance(design, str) or not design.strip().startswith("~"):
        report.error("design must be a nonempty R-style formula beginning with ~")
    elif "condition" not in design:
        report.error("design must contain the condition term used by the contrast")
    if numerator == denominator:
        report.error("contrast numerator and denominator must differ")
    for label, value in (("numerator", numerator), ("denominator", denominator)):
        if value not in conditions:
            report.error(
                f"contrast {label} {value!r} is not present in sample conditions "
                f"{sorted(conditions)}"
            )
    if reference != denominator:
        report.error("reference_level must equal contrast_denominator")
    if not isinstance(alpha, (int, float)) or isinstance(alpha, bool):
        report.error("adjusted_pvalue_threshold must be numeric")
    elif not 0 < float(alpha) <= 1:
        report.error("adjusted_pvalue_threshold must be greater than 0 and at most 1")
    if not isinstance(lfc_threshold, (int, float)) or isinstance(
        lfc_threshold, bool
    ):
        report.error("absolute_log2_fold_change_threshold must be numeric")
    elif float(lfc_threshold) < 0:
        report.error("absolute_log2_fold_change_threshold must be nonnegative")
    if shrinkage not in {"apeglm", "ashr", "none"}:
        report.error("lfc_shrinkage must be apeglm, ashr, or none")

    nominal = nested(
        data, ("downstream", "exploratory_nominal_analysis_allowed"), report
    )
    if not isinstance(nominal, bool):
        report.error("exploratory_nominal_analysis_allowed must be boolean")
    elif nominal:
        report.warn(
            "exploratory nominal analysis is enabled; primary and nominal "
            "outputs require explicit separation"
        )

    if not report.errors:
        report.passed("parameter document satisfies the agent contract checks")
    return data


def validation_has_pass(items: Any) -> bool:
    if not isinstance(items, list):
        return False
    for item in items:
        if isinstance(item, str) and "PASS" in item.upper():
            return True
        if isinstance(item, dict) and str(item.get("status", "")).upper() == "PASS":
            return True
    return False


def validate_status(
    path: Path,
    parameters: dict[str, Any] | None,
    report: Report,
) -> None:
    data = load_json(path, "analysis status", report)
    if not isinstance(data, dict):
        if data is not None:
            report.error("analysis status root must be a JSON object")
        return

    for location in find_placeholders(data):
        report.error(f"unresolved status placeholder at {location}")

    if data.get("schema_version") != 1:
        report.error("status.schema_version must equal 1")
    if data.get("mode") not in VALID_MODES:
        report.error(f"status.mode must be one of {sorted(VALID_MODES)}")
    if data.get("scope") != "conventional_bulk_paired_end_ATAC_seq":
        report.error("status.scope is outside this agent contract")
    if data.get("current_stage") not in STAGE_IDS:
        report.error("status.current_stage is not a recognized stage ID")

    if parameters is not None:
        if data.get("analysis_id") != parameters.get("analysis_id"):
            report.error("status and parameter analysis_id values do not match")

    decisions = data.get("human_decisions")
    if not isinstance(decisions, dict):
        report.error("status.human_decisions must be an object")
        decisions = {}
    missing_decisions = [key for key in DECISION_IDS if key not in decisions]
    if missing_decisions:
        report.error(
            "status is missing human decisions: " + ", ".join(missing_decisions)
        )

    approved_decisions: set[str] = set()
    for decision_id in DECISION_IDS:
        decision = decisions.get(decision_id)
        if not isinstance(decision, dict):
            continue
        decision_status = decision.get("status")
        if decision_status not in {"pending", "approved"}:
            report.error(
                f"human decision {decision_id} status must be pending or approved"
            )
        elif decision_status == "approved":
            approved_decisions.add(decision_id)
            for field in ("approved_by", "approved_at", "decision"):
                if not decision.get(field):
                    report.error(
                        f"approved human decision {decision_id} lacks {field}"
                    )
        else:
            report.warn(f"human decision remains pending: {decision_id}")

    stages = data.get("stages")
    if not isinstance(stages, list):
        report.error("status.stages must be an array")
        return
    observed_ids = [stage.get("id") for stage in stages if isinstance(stage, dict)]
    if observed_ids != STAGE_IDS:
        report.error(
            "status stage order differs from the workflow contract: "
            f"expected {STAGE_IDS}; observed {observed_ids}"
        )
        return

    statuses: dict[str, str] = {}
    for stage in stages:
        stage_id = stage["id"]
        stage_status = stage.get("status")
        statuses[stage_id] = str(stage_status)
        if stage_status not in VALID_STAGE_STATUSES:
            report.error(
                f"stage {stage_id} has invalid status {stage_status!r}"
            )
            continue
        if stage_status == "complete" and not validation_has_pass(
            stage.get("validation")
        ):
            report.error(
                f"complete stage {stage_id} has no recorded PASS validation"
            )

    active_statuses = {"approved", "running", "complete"}
    completed_statuses = {"complete", "not_applicable"}
    for stage_id, stage_status in statuses.items():
        if stage_status not in active_statuses:
            continue
        for dependency in DEPENDENCIES[stage_id]:
            if statuses.get(dependency) not in completed_statuses:
                report.error(
                    f"stage {stage_id} is {stage_status} but dependency "
                    f"{dependency} is {statuses.get(dependency)}"
                )
        for decision_id in STAGE_DECISIONS[stage_id]:
            if decision_id not in approved_decisions:
                report.error(
                    f"stage {stage_id} is {stage_status} without approved "
                    f"decision {decision_id}"
                )

    if data.get("mode") != "new_study_execution":
        running = [key for key, value in statuses.items() if value == "running"]
        if running:
            report.error(
                f"non-execution mode cannot contain running stages: {running}"
            )

    if not report.errors:
        report.passed("analysis status is structurally consistent")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Validate an AI-assisted bulk paired-end ATAC-seq configuration "
            "without modifying any files."
        )
    )
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--parameters", required=True, type=Path)
    parser.add_argument("--status", type=Path)
    parser.add_argument("--contract", type=Path, default=DEFAULT_CONTRACT)
    parser.add_argument(
        "--check-files",
        action="store_true",
        help="Require every recorded FASTQ path to exist and be nonempty.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Return a nonzero exit status when warnings are present.",
    )
    parser.add_argument(
        "--format",
        choices=("text", "json"),
        default="text",
        help="Select human-readable or machine-readable output.",
    )
    return parser.parse_args()


def emit(report: Report, output_format: str) -> None:
    if report.errors:
        overall_status = "FAIL"
    elif report.warnings:
        overall_status = "PASS_WITH_WARNINGS"
    else:
        overall_status = "PASS"
    payload = {
        "status": overall_status,
        "errors": report.errors,
        "warnings": report.warnings,
        "passes": report.passes,
        "counts": {
            "errors": len(report.errors),
            "warnings": len(report.warnings),
            "passes": len(report.passes),
        },
    }
    if output_format == "json":
        print(json.dumps(payload, indent=2, sort_keys=True))
        return

    print(f"status={payload['status']}")
    labels = {"errors": "ERROR", "warnings": "WARNING", "passes": "PASS"}
    for category in ("errors", "warnings", "passes"):
        for message in payload[category]:
            print(f"{labels[category]}: {message}")
    print(
        "summary="
        f"{len(report.errors)} errors, "
        f"{len(report.warnings)} warnings, "
        f"{len(report.passes)} passed checks"
    )


def main() -> int:
    args = parse_args()
    report = Report()

    validate_contract(args.contract, report)
    rows, conditions = validate_samples(args.samples, report, args.check_files)
    parameters = validate_parameters(
        args.parameters,
        rows,
        conditions,
        report,
    )
    if args.status is not None:
        validate_status(args.status, parameters, report)

    emit(report, args.format)
    if report.errors:
        return 1
    if args.strict and report.warnings:
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
