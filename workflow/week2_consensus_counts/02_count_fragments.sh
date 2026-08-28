#!/usr/bin/env bash
#
# Count paired ATAC-seq fragments in the consensus SAF with featureCounts.
# Refactored reference version: not the script that produced attempt1.

set -euo pipefail
export LC_ALL=C
umask 027

: "${PROJECT_ROOT:?Set PROJECT_ROOT to the repository root}"
: "${INPUT_ROOT:?Set INPUT_ROOT to the directory containing BAM and broadPeak inputs}"
TOOL_DIR="${TOOL_DIR:-}"
EXPECTED_HOST="${EXPECTED_HOST:-}"
ATTEMPT="${ATTEMPT:-refactored}"
MANIFEST="${MANIFEST:-${PROJECT_ROOT}/config/week2_samples.tsv}"
THREADS=2
EXPECTED_FEATURES="${EXPECTED_FEATURES:-112759}"

[[ "$ATTEMPT" =~ ^[A-Za-z0-9._-]+$ ]] || {
    printf 'ERROR: unsafe ATTEMPT name: %s\n' "$ATTEMPT" >&2
    exit 1
}

ATTEMPT_ROOT="${PROJECT_ROOT}/results/week2_consensus_counts/${ATTEMPT}"
SAF="${ATTEMPT_ROOT}/annotation/mature_memory_consensus.d0.saf"
COUNTS="${ATTEMPT_ROOT}/counts"
RAW_DIR="${COUNTS}/raw"
METADATA_DIR="${COUNTS}/metadata"
VALIDATION_DIR="${COUNTS}/validation"
RUNTIME="${COUNTS}/featurecounts_runtime"
LOG="${PROJECT_ROOT}/logs/week2_consensus_counts/${ATTEMPT}/featurecounts"
RAW="${RAW_DIR}/mature_memory_featurecounts.raw.tsv"
SUMMARY="${RAW}.summary"
ORDERED_SAMPLES="${METADATA_DIR}/ordered_samples.tsv"
VALIDATION="${VALIDATION_DIR}/featurecounts_validation.tsv"

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

if [[ -n "$EXPECTED_HOST" ]]; then
    [[ "$(hostname -s)" == "$EXPECTED_HOST" ]] ||
        die "expected host $EXPECTED_HOST; observed $(hostname -s)"
fi
[[ "$(pwd -P)" == "$PROJECT_ROOT" ]] ||
    die "run from exactly $PROJECT_ROOT"
[[ -s "$MANIFEST" ]] || die "missing sample manifest: $MANIFEST"
[[ -s "$SAF" ]] || die "missing consensus SAF: $SAF"
[[ ! -e "$COUNTS" && ! -L "$COUNTS" ]] ||
    die "refusing to overwrite counts path: $COUNTS"
[[ ! -e "$LOG" && ! -L "$LOG" ]] ||
    die "refusing to overwrite log path: $LOG"
[[ "$THREADS" -eq 2 ]] || die "featureCounts must use exactly two threads"

if [[ -n "$TOOL_DIR" ]]; then
    export PATH="${TOOL_DIR}:${PATH}"
fi
for command_name in awk featureCounts hostname ln mv pgrep readlink samtools sha256sum; do
    command -v "$command_name" >/dev/null ||
        die "required command not found: $command_name"
done

featurecounts_version="$(featureCounts -v 2>&1 | tr '\n' ' ')"
[[ "$featurecounts_version" =~ (^|[^0-9])v?2\.1\.1([^0-9.]|$) ]] ||
    die "featureCounts 2.1.1 is required; observed: $featurecounts_version"
samtools_version_output="$(samtools --version 2>&1)"
samtools_version="${samtools_version_output%%$'\n'*}"
[[ "$samtools_version" =~ (^|[[:space:]])1\.24($|[[:space:]]) ]] ||
    die "SAMtools 1.24 is required; observed: $samtools_version"

EXPECTED_HEADER=$'order\tsample_id\tdonor\tcell_type\tpairing_status\tbroadpeak_file\tbam_file'
IFS= read -r observed_header < "$MANIFEST"
[[ "$observed_header" == "$EXPECTED_HEADER" ]] ||
    die "unexpected sample-manifest header"
awk -F '\t' '
    NR == 1 { next }
    NF != 7 || $1 != NR - 1 || $7 ~ /\// || seen[$2]++ { exit 1 }
    $4 == "Mature_NK" { mature++ }
    $4 == "Memory_NK" { memory++ }
    END { if (NR != 11 || mature != 4 || memory != 6) exit 1 }
' "$MANIFEST" || die "sample manifest does not match the locked ten-sample design"

sample_ids=()
donors=()
cell_types=()
pairing_statuses=()
bam_paths=()
bam_names=()

while IFS=$'\t' read -r order sample_id donor cell_type pairing broadpeak bam; do
    bam_path="${INPUT_ROOT}/${bam}"
    [[ -s "$bam_path" ]] || die "missing BAM: $bam_path"
    [[ -s "${bam_path}.bai" ]] || die "missing BAI: ${bam_path}.bai"
    samtools quickcheck -v "$bam_path" ||
        die "samtools quickcheck failed: $bam_path"
    samtools view -H "$bam_path" |
        awk -F '\t' '$1 == "@HD" { for (i=2; i<=NF; i++) if ($i == "SO:coordinate") ok=1 } END { exit !ok }' ||
        die "BAM is not declared coordinate sorted: $bam_path"

    sample_ids+=("$sample_id")
    donors+=("$donor")
    cell_types+=("$cell_type")
    pairing_statuses+=("$pairing")
    bam_paths+=("$bam_path")
    bam_names+=("$bam")
done < <(tail -n +2 "$MANIFEST")

[[ ${#bam_paths[@]} -eq 10 ]] || die "expected exactly ten BAM inputs"

for process_name in featureCounts bedtools R Rscript macs3 bowtie2 cutadapt fastqc multiqc picard; do
    pgrep -u "$(id -u)" -x "$process_name" >/dev/null 2>&1 &&
        die "conflicting analytical process is running: $process_name"
done

if ! awk -F '\t' -v expected="$EXPECTED_FEATURES" '
    NR == 1 {
        if ($0 != "GeneID\tChr\tStart\tEnd\tStrand") exit 1
        next
    }
    NF != 5 || $1 != sprintf("peak_%06d", NR - 1) ||
    $3 !~ /^[0-9]+$/ || $4 !~ /^[0-9]+$/ ||
    $3 < 1 || $4 < $3 || $5 != "." { exit 1 }
    END { if (NR - 1 != expected) exit 1 }
' "$SAF"; then
    die "SAF structure or feature count is invalid"
fi

mkdir -p "$RAW_DIR" "$METADATA_DIR" "$VALIDATION_DIR" "$RUNTIME" "$LOG"
RUN_LOG="${LOG}/02_count_fragments_refactored.log"
exec > >(tee "$RUN_LOG") 2>&1

printf 'featureCounts=%s\n' "$featurecounts_version"
printf 'samtools=%s\n' "$samtools_version"
printf 'threads=%s\n' "$THREADS"

{
    printf 'order\tsample_id\tdonor\tcell_type\tpairing_status\tbam_path\n'
    for index in "${!bam_paths[@]}"; do
        printf '%d\t%s\t%s\t%s\t%s\t%s\n' \
            "$((index + 1))" \
            "${sample_ids[$index]}" \
            "${donors[$index]}" \
            "${cell_types[$index]}" \
            "${pairing_statuses[$index]}" \
            "${bam_paths[$index]}"
    done
} > "$ORDERED_SAMPLES"

ln -s "$SAF" "${RUNTIME}/mature_memory_consensus.d0.saf"
runtime_bams=()
for index in "${!bam_paths[@]}"; do
    ln -s "${bam_paths[$index]}" "${RUNTIME}/${bam_names[$index]}"
    ln -s "${bam_paths[$index]}.bai" "${RUNTIME}/${bam_names[$index]}.bai"
    [[ "$(readlink -f "${RUNTIME}/${bam_names[$index]}")" == "${bam_paths[$index]}" ]] ||
        die "runtime BAM link target mismatch: ${bam_names[$index]}"
    runtime_bams+=("${bam_names[$index]}")
done

printf 'command=featureCounts -F SAF -a mature_memory_consensus.d0.saf '
printf -- '-o mature_memory_featurecounts.raw.tsv -p --countReadPairs '
printf -- '-B -C -s 0 -T 2 %s\n' "${runtime_bams[*]}"

(
    cd "$RUNTIME"
    featureCounts \
        -F SAF \
        -a mature_memory_consensus.d0.saf \
        -o mature_memory_featurecounts.raw.tsv \
        -p \
        --countReadPairs \
        -B \
        -C \
        -s 0 \
        -T "$THREADS" \
        "${runtime_bams[@]}"
)

RUNTIME_RAW="${RUNTIME}/mature_memory_featurecounts.raw.tsv"
RUNTIME_SUMMARY="${RUNTIME_RAW}.summary"
[[ -s "$RUNTIME_RAW" && -s "$RUNTIME_SUMMARY" ]] ||
    die "featureCounts did not create both raw and summary outputs"

expected_bams_csv="$(IFS=,; printf '%s' "${bam_names[*]}")"
if ! awk -F '\t' -v rows="$EXPECTED_FEATURES" -v names="$expected_bams_csv" '
    function base(path, n, parts) { n=split(path, parts, "/"); return parts[n] }
    BEGIN { ok=1; split(names, expected, ",") }
    /^#/ { comments++; next }
    !header {
        header=1
        if (NF != 16 || $1 != "Geneid" || $2 != "Chr" ||
            $3 != "Start" || $4 != "End" || $5 != "Strand" ||
            $6 != "Length") ok=0
        for (i=1; i<=10; i++) if (base($(i+6)) != expected[i]) ok=0
        next
    }
    {
        data++
        if (NF != 16 || seen[$1]++) ok=0
        for (i=7; i<=16; i++) if ($i !~ /^[0-9]+$/) ok=0
    }
    END { if (comments < 1 || data != rows || !ok) exit 1 }
' "$RUNTIME_RAW"; then
    die "raw featureCounts output failed structure, order, or count validation"
fi

if ! awk -F '\t' -v expected="$EXPECTED_FEATURES" '
    NR == FNR {
        if (FNR == 1) next
        saf_rows++
        id[saf_rows]=$1; chr[saf_rows]=$2; start[saf_rows]=$3
        end[saf_rows]=$4; strand[saf_rows]=$5
        next
    }
    /^#/ { next }
    !header { header=1; next }
    {
        rows++
        if ($1 != id[rows] || $2 != chr[rows] || $3 != start[rows] ||
            $4 != end[rows] || $5 != strand[rows] ||
            $6 != $4 - $3 + 1) mismatch=1
    }
    END {
        if (rows != expected || saf_rows != expected || mismatch) exit 1
    }
' "$SAF" "$RUNTIME_RAW"; then
    die "raw feature annotation does not exactly match the consensus SAF"
fi

if ! awk -F '\t' -v names="$expected_bams_csv" '
    function base(path, n, parts) { n=split(path, parts, "/"); return parts[n] }
    BEGIN { ok=1; split(names, expected, ",") }
    NR == 1 {
        if (NF != 11 || $1 != "Status") ok=0
        for (i=1; i<=10; i++) if (base($(i+1)) != expected[i]) ok=0
        next
    }
    {
        if (NF != 11 || seen[$1]++) ok=0
        if ($1 == "Assigned") assigned++
        for (i=2; i<=11; i++) if ($i !~ /^[0-9]+$/) ok=0
    }
    END { if (assigned != 1 || !ok) exit 1 }
' "$RUNTIME_SUMMARY"; then
    die "featureCounts summary failed structure or sample-order validation"
fi

if ! awk -F '\t' '
    NR == FNR {
        if ($0 ~ /^#/) next
        if (!raw_header++) next
        for (i=1; i<=10; i++) matrix[i] += $(i+6)
        next
    }
    FNR == 1 { next }
    $1 == "Assigned" {
        assigned++
        for (i=1; i<=10; i++) if (matrix[i] != $(i+1)) mismatch=1
    }
    END { if (assigned != 1 || mismatch) exit 1 }
' "$RUNTIME_RAW" "$RUNTIME_SUMMARY"; then
    die "matrix column sums do not match the summary Assigned row"
fi

mv "$RUNTIME_RAW" "$RAW"
mv "$RUNTIME_SUMMARY" "$SUMMARY"

{
    printf 'check\tvalue\tstatus\n'
    printf 'sample_count\t10\tPASS\n'
    printf 'feature_rows\t%s\tPASS\n' "$EXPECTED_FEATURES"
    printf 'sample_order_matches\tYES\tPASS\n'
    printf 'counts_are_nonnegative_integers\tYES\tPASS\n'
    printf 'raw_annotation_matches_SAF\tYES\tPASS\n'
    printf 'matrix_sums_match_Assigned\tYES\tPASS\n'
    printf 'featureCounts_options\t-F SAF -p --countReadPairs -B -C -s 0 -T 2\tPASS\n'
} > "$VALIDATION"

sha256sum "$ORDERED_SAMPLES" "$RAW" "$SUMMARY" "$VALIDATION" \
    > "${LOG}/02_count_fragments_refactored.sha256"

printf 'raw_output=%s\n' "$RAW"
printf 'summary=%s\n' "$SUMMARY"
printf 'clean_matrix_created=NO\n'
printf 'DESeq2_executed=NO\n'
printf 'status=PASS\n'
