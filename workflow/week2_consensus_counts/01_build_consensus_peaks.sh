#!/usr/bin/env bash
#
# Build the all-inclusive Mature/Memory consensus peak universe and SAF file.
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

[[ "$ATTEMPT" =~ ^[A-Za-z0-9._-]+$ ]] || {
    printf 'ERROR: unsafe ATTEMPT name: %s\n' "$ATTEMPT" >&2
    exit 1
}

OUT="${PROJECT_ROOT}/results/week2_consensus_counts/${ATTEMPT}"
LOG="${PROJECT_ROOT}/logs/week2_consensus_counts/${ATTEMPT}/consensus"
INTERMEDIATE="${OUT}/intermediate"
CONSENSUS="${OUT}/consensus"
ANNOTATION="${OUT}/annotation"
MATURE_BED="${CONSENSUS}/mature_nk_union.d0.bed"
MEMORY_BED="${CONSENSUS}/memory_nk_union.d0.bed"
FINAL_BED="${CONSENSUS}/mature_memory_consensus.d0.bed"
FINAL_SAF="${ANNOTATION}/mature_memory_consensus.d0.saf"

die() {
    printf 'ERROR: %s\n' "$*" >&2
    exit 1
}

rows() {
    awk 'END { print NR + 0 }' "$1"
}

validate_bed3() {
    local file=$1
    [[ -s "$file" ]] || die "empty BED file: $file"
    awk '
        NF != 3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ ||
        $2 < 0 || $3 <= $2 || $1 ~ /^chr/ { exit 1 }
    ' "$file" || die "invalid BED3 structure or coordinates: $file"
    sort -c -k1,1 -k2,2n -k3,3n "$file" ||
        die "BED3 is not coordinate sorted: $file"
}

validate_merged() {
    local file=$1
    awk '
        NR > 1 && $1 == previous_chr && $2 <= previous_end { exit 1 }
        { previous_chr = $1; previous_end = $3 }
    ' "$file" || die "overlapping or book-ended intervals remain: $file"
}

if [[ -n "$EXPECTED_HOST" ]]; then
    [[ "$(hostname -s)" == "$EXPECTED_HOST" ]] ||
        die "expected host $EXPECTED_HOST; observed $(hostname -s)"
fi
[[ "$(pwd -P)" == "$PROJECT_ROOT" ]] ||
    die "run from exactly $PROJECT_ROOT"
[[ -s "$MANIFEST" ]] || die "missing sample manifest: $MANIFEST"
[[ ! -e "$OUT" && ! -L "$OUT" ]] ||
    die "refusing to overwrite output path: $OUT"
[[ ! -e "$LOG" && ! -L "$LOG" ]] ||
    die "refusing to overwrite log path: $LOG"

if [[ -n "$TOOL_DIR" ]]; then
    export PATH="${TOOL_DIR}:${PATH}"
fi
for command_name in awk bedtools date hostname sha256sum sort; do
    command -v "$command_name" >/dev/null ||
        die "required command not found: $command_name"
done

EXPECTED_HEADER=$'order\tsample_id\tdonor\tcell_type\tpairing_status\tbroadpeak_file\tbam_file'
IFS= read -r observed_header < "$MANIFEST"
[[ "$observed_header" == "$EXPECTED_HEADER" ]] ||
    die "unexpected sample-manifest header"

awk -F '\t' '
    NR == 1 { next }
    NF != 7 || $1 != NR - 1 || $2 !~ /^SRR[0-9]+$/ ||
    ($4 != "Mature_NK" && $4 != "Memory_NK") ||
    $6 ~ /\// || $7 ~ /\// || seen[$2]++ { exit 1 }
    $4 == "Mature_NK" { mature++ }
    $4 == "Memory_NK" { memory++ }
    END { if (NR != 11 || mature != 4 || memory != 6) exit 1 }
' "$MANIFEST" || die "sample manifest must contain the locked 4 Mature + 6 Memory order"

mkdir -p "$INTERMEDIATE" "$CONSENSUS" "$ANNOTATION" "$LOG"
RUN_LOG="${LOG}/01_build_consensus_peaks_refactored.log"
exec > >(tee "$RUN_LOG") 2>&1

printf 'bedtools=%s\n' "$(bedtools --version)"
printf 'attempt=%s\n' "$ATTEMPT"

mature_inputs=()
memory_inputs=()

while IFS=$'\t' read -r order sample_id donor cell_type pairing broadpeak bam; do
    source_file="${INPUT_ROOT}/${broadpeak}"
    extracted="${INTERMEDIATE}/${sample_id}.broadPeak.bed3"
    sorted="${INTERMEDIATE}/${sample_id}.broadPeak.sorted.bed3"

    [[ -s "$source_file" ]] || die "missing broadPeak input: $source_file"
    awk '
        NF != 9 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ ||
        $2 < 0 || $3 <= $2 || $1 ~ /^chr/ { exit 1 }
    ' "$source_file" || die "invalid broadPeak input: $source_file"

    awk 'BEGIN { OFS = "\t" } { print $1, $2, $3 }' \
        "$source_file" > "$extracted"
    sort -k1,1 -k2,2n -k3,3n "$extracted" > "$sorted"
    validate_bed3 "$sorted"
    [[ "$(rows "$source_file")" == "$(rows "$sorted")" ]] ||
        die "peak rows changed during extraction: $sample_id"

    if [[ "$cell_type" == "Mature_NK" ]]; then
        mature_inputs+=("$sorted")
    else
        memory_inputs+=("$sorted")
    fi
done < <(tail -n +2 "$MANIFEST")

[[ ${#mature_inputs[@]} -eq 4 && ${#memory_inputs[@]} -eq 6 ]] ||
    die "internal sample-group count mismatch"

merge_condition() {
    local output=$1
    shift
    local combined="${output%.bed}.all.sorted.bed"
    awk 'BEGIN { OFS = "\t" } { print $1, $2, $3 }' "$@" |
        sort -k1,1 -k2,2n -k3,3n > "$combined"
    bedtools merge -d 0 -i "$combined" > "$output"
    validate_bed3 "$output"
    validate_merged "$output"
}

merge_condition "$MATURE_BED" "${mature_inputs[@]}"
merge_condition "$MEMORY_BED" "${memory_inputs[@]}"

condition_unions="${INTERMEDIATE}/mature_memory.condition_unions.sorted.bed3"
awk 'BEGIN { OFS = "\t" } { print $1, $2, $3 }' \
    "$MATURE_BED" "$MEMORY_BED" |
    sort -k1,1 -k2,2n -k3,3n > "$condition_unions"
bedtools merge -d 0 -i "$condition_unions" > "$FINAL_BED"
validate_bed3 "$FINAL_BED"
validate_merged "$FINAL_BED"

awk '
    BEGIN { OFS = "\t"; print "GeneID", "Chr", "Start", "End", "Strand" }
    { printf "peak_%06d\t%s\t%d\t%d\t.\n", NR, $1, $2 + 1, $3 }
' "$FINAL_BED" > "$FINAL_SAF"

awk -F '\t' '
    NR == FNR {
        if (FNR == 1) {
            if ($0 != "GeneID\tChr\tStart\tEnd\tStrand") exit 1
            next
        }
        n++
        id[n] = $1; chr[n] = $2; start[n] = $3; end[n] = $4
        next
    }
    {
        row = FNR
        expected = sprintf("peak_%06d", row)
        if (id[row] != expected || chr[row] != $1 ||
            start[row] != $2 + 1 || end[row] != $3 ||
            end[row] - start[row] + 1 != $3 - $2) exit 1
    }
    END { if (n != FNR) exit 1 }
' "$FINAL_SAF" "$FINAL_BED" || die "BED-to-SAF equivalence check failed"

sha256sum "$MATURE_BED" "$MEMORY_BED" "$FINAL_BED" "$FINAL_SAF" \
    > "${LOG}/01_build_consensus_peaks_refactored.sha256"

printf 'mature_union_regions=%s\n' "$(rows "$MATURE_BED")"
printf 'memory_union_regions=%s\n' "$(rows "$MEMORY_BED")"
printf 'final_consensus_regions=%s\n' "$(rows "$FINAL_BED")"
printf 'status=PASS\n'
