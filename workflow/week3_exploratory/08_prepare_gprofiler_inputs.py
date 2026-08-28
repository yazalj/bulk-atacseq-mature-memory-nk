#!/usr/bin/env python3
"""Prepare reproducible custom-background g:Profiler inputs for Week 3.

This script does not rerun DESeq2, apeglm, or HOMER. It maps the protected
59,186-peak universe and the approved nominal candidate sets to the nearest
GENCODE v50 transcript TSS, then prepares direction-aware g:Profiler queries.
Nearest-gene assignments are hypotheses, not validated regulatory links.
"""

from __future__ import annotations

import bisect
import csv
import gzip
import hashlib
import json
import os
import sys
from collections import defaultdict
from pathlib import Path


EXPECTED_SHA256 = {
    "primary_unshrunken_results": "d02676374a1b58c1913254e8c9b3f0c1dec4abfdc2ecf31756742871e1f01893",
    "primary_candidate_table": "9218ff76ccf340fa027faa1ddd9ec09969b0f86035be704cdb2838cb14b8092e",
    "gencode_v50_gtf": "89bbad69a8c89fee5fadec0a2f14a098752d27a4ef7f24e4de9bac681e1b18f4",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_gtf_attributes(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for item in text.rstrip("; ").split(";"):
        item = item.strip()
        if not item:
            continue
        key, _, value = item.partition(" ")
        out[key] = value.strip().strip('"')
    return out


def gencode_seqname(chrom: str) -> str:
    if chrom in {"M", "MT", "chrM"}:
        return "chrM"
    return chrom if chrom.startswith("chr") else f"chr{chrom}"


def truthy(value: str) -> bool:
    return value.strip().upper() == "TRUE"


def write_tsv(path: Path, fieldnames: list[str], rows) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_list(path: Path, values: list[str]) -> None:
    path.write_text("\n".join(values) + "\n", encoding="utf-8")


def main() -> int:
    if len(sys.argv) != 3:
        raise SystemExit("Usage: 14_prepare_gprofiler_background_enrichment_attempt2.py <project_root> <gencode_v50_gtf_gz>")

    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["OPENBLAS_NUM_THREADS"] = "1"
    os.environ["MKL_NUM_THREADS"] = "1"

    project = Path(sys.argv[1]).resolve(strict=True)
    gtf_path = Path(sys.argv[2]).resolve(strict=True)
    result_root = project / "results" / "week3_downstream" / "gprofiler" / "attempt2"
    if result_root.exists():
        raise SystemExit(f"Refusing to overwrite existing output: {result_root}")

    table_dir = result_root / "tables"
    gene_dir = result_root / "gene_lists"
    provenance_dir = result_root / "provenance"
    for directory in (table_dir, gene_dir, provenance_dir):
        directory.mkdir(parents=True, exist_ok=False)

    ordinary_path = project / "results" / "week2_differential_accessibility" / "attempt3" / "tables" / "deseq2_results_unshrunken.tsv"
    candidate_path = project / "results" / "week2_differential_accessibility" / "nominal_sensitivity" / "attempt1" / "tables" / "raw_pvalue_candidates_minimum_samples_5.tsv"
    inputs = {
        "primary_unshrunken_results": ordinary_path,
        "primary_candidate_table": candidate_path,
        "gencode_v50_gtf": gtf_path,
    }
    for label, path in inputs.items():
        if not path.is_file():
            raise SystemExit(f"Missing required input: {label}: {path}")

    hash_rows = []
    for label, path in inputs.items():
        actual = sha256(path)
        expected = EXPECTED_SHA256[label]
        hash_rows.append({"input_id": label, "path": str(path), "expected_sha256": expected, "actual_sha256": actual, "match": str(actual == expected).upper()})
        if actual != expected:
            raise SystemExit(f"Protected/reference hash mismatch: {label}")
    write_tsv(provenance_dir / "input_sha256_verification.tsv", list(hash_rows[0]), hash_rows)

    transcripts_by_chr: dict[str, dict[int, list[dict[str, object]]]] = defaultdict(lambda: defaultdict(list))
    with gzip.open(gtf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            attrs = parse_gtf_attributes(fields[8])
            strand = fields[6]
            start = int(fields[3])
            end = int(fields[4])
            tss = start if strand == "+" else end
            gene_id = attrs.get("gene_id", "")
            transcript_id = attrs.get("transcript_id", "")
            if not gene_id or not transcript_id:
                continue
            transcripts_by_chr[fields[0]][tss].append(
                {
                    "tss": tss,
                    "strand": strand,
                    "gene_id": gene_id,
                    "gene_id_unversioned": gene_id.split(".", 1)[0],
                    "gene_symbol": attrs.get("gene_name", ""),
                    "gene_type": attrs.get("gene_type", ""),
                    "transcript_id": transcript_id,
                    "transcript_type": attrs.get("transcript_type", ""),
                }
            )

    positions_by_chr = {chrom: sorted(groups) for chrom, groups in transcripts_by_chr.items()}
    if not positions_by_chr:
        raise SystemExit("No GENCODE transcript TSS records were parsed")

    def nearest(chrom: str, start: int, end: int) -> dict[str, object]:
        seqname = gencode_seqname(chrom)
        positions = positions_by_chr.get(seqname)
        if not positions:
            raise ValueError(f"No GENCODE transcript TSS records for {seqname}")
        left = bisect.bisect_left(positions, start)
        right = bisect.bisect_right(positions, end)
        candidate_positions: list[int]
        if right > left:
            candidate_positions = positions[left:right]
            min_distance = 0
        else:
            candidate_positions = []
            if left > 0:
                candidate_positions.append(positions[left - 1])
            if left < len(positions):
                candidate_positions.append(positions[left])
            distances = [start - pos if pos < start else pos - end for pos in candidate_positions]
            min_distance = min(distances)
            candidate_positions = [pos for pos, distance in zip(candidate_positions, distances) if distance == min_distance]
        tied = [record for pos in candidate_positions for record in transcripts_by_chr[seqname][pos]]
        tied.sort(
            key=lambda record: (
                record["gene_type"] != "protein_coding",
                not bool(record["gene_symbol"]),
                str(record["gene_symbol"]),
                str(record["transcript_id"]),
            )
        )
        chosen = dict(tied[0])
        tss = int(chosen["tss"])
        strand = str(chosen["strand"])
        if start <= tss <= end:
            signed_distance = 0
        elif strand == "+":
            signed_distance = end - tss if end < tss else start - tss
        else:
            signed_distance = tss - end if end < tss else tss - start
        chosen.update(
            {
                "nearest_tss_distance_bp": min_distance,
                "signed_distance_to_tss_bp": signed_distance,
                "nearest_tss_tie_transcripts": len(tied),
                "nearest_tss_tie_genes": len({str(record["gene_id"]) for record in tied}),
            }
        )
        return chosen

    with ordinary_path.open("r", newline="", encoding="utf-8") as handle:
        ordinary_rows = list(csv.DictReader(handle, delimiter="\t"))
    with candidate_path.open("r", newline="", encoding="utf-8") as handle:
        candidate_rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(ordinary_rows) != 59186 or len(candidate_rows) != 3426:
        raise SystemExit(f"Unexpected protected row counts: primary={len(ordinary_rows)}, candidates={len(candidate_rows)}")
    if len({row["peak_id"] for row in ordinary_rows}) != len(ordinary_rows):
        raise SystemExit("Primary result table contains duplicate peak IDs")

    candidate_by_id = {row["peak_id"]: row for row in candidate_rows}
    broader_ids = {row["peak_id"] for row in candidate_rows if truthy(row["candidate_raw_p_and_ordinary_lfc"])}
    strict_ids = {row["peak_id"] for row in candidate_rows if truthy(row["candidate_raw_p_and_apeglm_lfc"])}
    if len(broader_ids) != 1562 or len(strict_ids) != 21:
        raise SystemExit(f"Approved target sets did not reproduce: broader={len(broader_ids)}, strict={len(strict_ids)}")

    mapping_fields = [
        "peak_id", "Chr", "Start", "End", "log2FoldChange", "pvalue", "padj",
        "nearest_gene_id", "nearest_gene_id_unversioned", "nearest_gene_symbol", "nearest_gene_type",
        "nearest_transcript_id", "nearest_transcript_type", "nearest_transcript_strand", "nearest_tss",
        "nearest_tss_distance_bp", "signed_distance_to_tss_bp", "nearest_tss_tie_transcripts", "nearest_tss_tie_genes",
        "is_broader_ordinary_candidate", "is_strict_apeglm_candidate", "ordinary_direction",
    ]
    all_mappings: list[dict[str, object]] = []
    broader_mappings: list[dict[str, object]] = []
    strict_mappings: list[dict[str, object]] = []
    unmappable_rows: list[dict[str, object]] = []
    for row in ordinary_rows:
        peak_id = row["peak_id"]
        try:
            hit = nearest(row["Chr"], int(row["Start"]), int(row["End"]))
        except ValueError:
            unmappable_rows.append(
                {
                    "peak_id": peak_id,
                    "Chr": row["Chr"],
                    "Start": row["Start"],
                    "End": row["End"],
                    "is_broader_ordinary_candidate": str(peak_id in broader_ids).upper(),
                    "is_strict_apeglm_candidate": str(peak_id in strict_ids).upper(),
                    "exclusion_reason": "No transcript TSS for contig in GENCODE v50 primary-assembly annotation",
                }
            )
            continue
        fold_change = float(row["log2FoldChange"])
        mapped: dict[str, object] = {
            "peak_id": peak_id,
            "Chr": row["Chr"],
            "Start": row["Start"],
            "End": row["End"],
            "log2FoldChange": row["log2FoldChange"],
            "pvalue": row["pvalue"],
            "padj": row["padj"],
            "nearest_gene_id": hit["gene_id"],
            "nearest_gene_id_unversioned": hit["gene_id_unversioned"],
            "nearest_gene_symbol": hit["gene_symbol"],
            "nearest_gene_type": hit["gene_type"],
            "nearest_transcript_id": hit["transcript_id"],
            "nearest_transcript_type": hit["transcript_type"],
            "nearest_transcript_strand": hit["strand"],
            "nearest_tss": hit["tss"],
            "nearest_tss_distance_bp": hit["nearest_tss_distance_bp"],
            "signed_distance_to_tss_bp": hit["signed_distance_to_tss_bp"],
            "nearest_tss_tie_transcripts": hit["nearest_tss_tie_transcripts"],
            "nearest_tss_tie_genes": hit["nearest_tss_tie_genes"],
            "is_broader_ordinary_candidate": str(peak_id in broader_ids).upper(),
            "is_strict_apeglm_candidate": str(peak_id in strict_ids).upper(),
            "ordinary_direction": "Memory_NK_more_accessible" if fold_change > 0 else "Mature_NK_more_accessible",
        }
        all_mappings.append(mapped)
        if peak_id in broader_ids:
            broader_mappings.append(mapped)
        if peak_id in strict_ids:
            strict_mappings.append(mapped)

    if len(all_mappings) + len(unmappable_rows) != 59186 or len(strict_mappings) != 21:
        raise SystemExit("Nearest-gene mapping accounting or strict-set row counts failed")

    background_gz = table_dir / "primary_59186_peak_to_nearest_gencode_gene.tsv.gz"
    with gzip.open(background_gz, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=mapping_fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_mappings)
    write_tsv(table_dir / "broader_1562_peak_to_nearest_gencode_gene.tsv", mapping_fields, broader_mappings)
    write_tsv(table_dir / "strict_21_peak_to_nearest_gencode_gene.tsv", mapping_fields, strict_mappings)
    if unmappable_rows:
        write_tsv(table_dir / "annotation_unmappable_peaks.tsv", list(unmappable_rows[0]), unmappable_rows)

    background_genes = sorted({str(row["nearest_gene_id_unversioned"]) for row in all_mappings})
    broader_genes = sorted({str(row["nearest_gene_id_unversioned"]) for row in broader_mappings})
    memory_genes = sorted({str(row["nearest_gene_id_unversioned"]) for row in broader_mappings if row["ordinary_direction"] == "Memory_NK_more_accessible"})
    mature_genes = sorted({str(row["nearest_gene_id_unversioned"]) for row in broader_mappings if row["ordinary_direction"] == "Mature_NK_more_accessible"})
    strict_genes = sorted({str(row["nearest_gene_id_unversioned"]) for row in strict_mappings})
    if not set(broader_genes).issubset(background_genes) or not set(strict_genes).issubset(background_genes):
        raise SystemExit("One or more query genes are absent from the custom background")

    gene_sets = {
        "background_primary_59186_nearest_genes.txt": background_genes,
        "broader_1562_all_nearest_genes.txt": broader_genes,
        "broader_1562_memory_more_nearest_genes.txt": memory_genes,
        "broader_1562_mature_more_nearest_genes.txt": mature_genes,
        "strict_21_nearest_genes.txt": strict_genes,
    }
    for name, values in gene_sets.items():
        write_list(gene_dir / name, values)

    payload = {
        "organism": "hsapiens",
        "query": {
            "Memory_more_broader_ordinary": memory_genes,
            "Mature_more_broader_ordinary": mature_genes,
            "Strict_21_apeglm": strict_genes,
        },
        "sources": ["GO:BP", "GO:MF", "GO:CC", "REAC"],
        "user_threshold": 0.05,
        "all_results": True,
        "ordered": False,
        "combined": False,
        "measure_underrepresentation": False,
        "no_iea": False,
        "domain_scope": "custom",
        "background": background_genes,
        "significance_threshold_method": "fdr",
        "no_evidences": True,
    }
    (result_root / "gprofiler_request.json").write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    summary = [
        {"metric": "primary_peaks_input", "value": len(ordinary_rows)},
        {"metric": "primary_peaks_mapped", "value": len(all_mappings)},
        {"metric": "primary_peaks_annotation_unmappable", "value": len(unmappable_rows)},
        {"metric": "custom_background_unique_genes", "value": len(background_genes)},
        {"metric": "broader_candidate_peaks_input", "value": len(broader_ids)},
        {"metric": "broader_candidate_peaks_mapped", "value": len(broader_mappings)},
        {"metric": "broader_unique_genes", "value": len(broader_genes)},
        {"metric": "memory_more_broader_unique_genes", "value": len(memory_genes)},
        {"metric": "mature_more_broader_unique_genes", "value": len(mature_genes)},
        {"metric": "strict_candidate_peaks", "value": len(strict_mappings)},
        {"metric": "strict_unique_genes", "value": len(strict_genes)},
        {"metric": "nearest_gene_assignment_status", "value": "HYPOTHESIS_ONLY"},
    ]
    write_tsv(provenance_dir / "preparation_metrics.tsv", ["metric", "value"], summary)
    (result_root / "STATUS.txt").write_text(
        "STATUS=PASS\n"
        "PURPOSE=CUSTOM_BACKGROUND_GPROFILER_INPUT_PREPARATION_WITH_EXPLICIT_UNMAPPABLE_CONTIG_EXCLUSION\n"
        "THREADS=1\n"
        "CAVEAT=Nearest-gene assignments are hypotheses, not validated regulatory links.\n",
        encoding="utf-8",
    )
    print(json.dumps({row["metric"]: row["value"] for row in summary}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
