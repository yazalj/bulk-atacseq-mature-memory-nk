#!/usr/bin/env bash
set -euo pipefail

: "${PROJECT_ROOT:?Set PROJECT_ROOT using the current shell path syntax}"
: "${HOMER_ROOT:?Set HOMER_ROOT to the HOMER installation root}"
: "${GENOME_FASTA:?Set GENOME_FASTA to the uncompressed GRCh38 FASTA}"
: "${GENOME_FASTA_GZ:?Set GENOME_FASTA_GZ to its gzip-compressed source}"

project_root="$PROJECT_ROOT"
homer_root="$HOMER_ROOT"
fasta_gz="$GENOME_FASTA_GZ"
fasta="$GENOME_FASTA"
validated_input_root="${HOMER_INPUT_ROOT:-$project_root/results/week3_downstream/homer/input}"
result_root="${HOMER_OUTPUT_ROOT:-$project_root/results/week3_downstream/homer/public_reproduction}"
log_root="${HOMER_LOG_ROOT:-$project_root/logs/week3_downstream/homer/public_reproduction}"

for required_directory in "$result_root/input" "$result_root/preflight" "$result_root/homer_output" "$log_root"; do
    [[ -d "$required_directory" ]] || { printf 'Missing pre-created Windows directory: %s\n' "$required_directory" >&2; exit 1; }
done
if [[ -n "$(find "$result_root" -type f -print -quit 2>/dev/null)" || -n "$(find "$log_root" -type f -print -quit 2>/dev/null)" ]]; then
    printf 'Refusing to overwrite an existing HOMER Attempt 4 file.\n' >&2
    exit 1
fi

exec > >(tee "$log_root/execution.log") 2>&1
export PATH="$homer_root/bin:$PATH"
export LC_ALL=C

target_source="$validated_input_root/broader_1562_homer_sensitivity_candidates_homer_contigs.bed"
background_source="$validated_input_root/primary_59186_tested_peaks_background_homer_contigs.bed"
target_bed="$result_root/input/broader_1562_homer_sensitivity_candidates_homer_contigs.bed"
background_bed="$result_root/input/primary_59186_tested_peaks_background_homer_contigs.bed"

printf 'WEEK3_HOMER_RUN=public_reproduction\n'
printf 'START_UTC=%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
printf 'HOST=%s\n' "$(hostname)"
printf 'CYGWIN=%s\n' "$(uname -srvmo)"
printf 'HOMER_ROOT=%s\n' "$homer_root"
printf 'FASTA=%s\n' "$fasta"

for required in "$target_source" "$background_source" "$fasta_gz" "$fasta" "$homer_root/bin/findMotifsGenome.pl" "$homer_root/bin/mergePeaks.exe" "$homer_root/bin/homerTools.exe" "$homer_root/bin/homer2.exe"; do
    [[ -s "$required" ]] || { printf 'Missing non-empty required file: %s\n' "$required" >&2; exit 1; }
done
[[ "$(md5sum "$fasta_gz" | awk '{print $1}')" == 'da1a11258be075cfa7af718162c894e7' ]] || { printf 'GENCODE FASTA MD5 mismatch.\n' >&2; exit 1; }
gzip -t "$fasta_gz"
perl -c "$homer_root/bin/findMotifsGenome.pl"
perl -c "$homer_root/bin/HomerSVGLogo.pm"

cp "$target_source" "$target_bed"
cp "$background_source" "$background_bed"

target_n=$(wc -l < "$target_bed")
background_n=$(wc -l < "$background_bed")
[[ "$target_n" -eq 1562 && "$background_n" -eq 59186 ]] || { printf 'Unexpected target/background counts.\n' >&2; exit 1; }
awk '$2 < 0 || $3 <= $2 {bad++} END {exit bad > 0}' "$target_bed" || { printf 'Invalid target coordinates.\n' >&2; exit 1; }
awk '$2 < 0 || $3 <= $2 {bad++} END {exit bad > 0}' "$background_bed" || { printf 'Invalid background coordinates.\n' >&2; exit 1; }
[[ -z "$(cut -f4 "$target_bed" | sort | uniq -d)" ]] || { printf 'Duplicate target IDs.\n' >&2; exit 1; }
[[ -z "$(cut -f4 "$background_bed" | sort | uniq -d)" ]] || { printf 'Duplicate background IDs.\n' >&2; exit 1; }

cut -f4 "$target_bed" | sort -u > "$result_root/preflight/target_ids.txt"
cut -f4 "$background_bed" | sort -u > "$result_root/preflight/background_ids.txt"
comm -23 "$result_root/preflight/target_ids.txt" "$result_root/preflight/background_ids.txt" > "$result_root/preflight/target_ids_missing_from_background.txt"
[[ ! -s "$result_root/preflight/target_ids_missing_from_background.txt" ]] || { printf 'Target is not a subset of background.\n' >&2; exit 1; }

grep '^>' "$fasta" | sed 's/^>//; s/ .*//' | sort -u > "$result_root/preflight/fasta_contigs.txt"
cut -f1 "$target_bed" | sort -u > "$result_root/preflight/target_contigs.txt"
cut -f1 "$background_bed" | sort -u > "$result_root/preflight/background_contigs.txt"
comm -23 "$result_root/preflight/target_contigs.txt" "$result_root/preflight/fasta_contigs.txt" > "$result_root/preflight/target_contigs_missing_from_fasta.txt"
comm -23 "$result_root/preflight/background_contigs.txt" "$result_root/preflight/fasta_contigs.txt" > "$result_root/preflight/background_contigs_missing_from_fasta.txt"
[[ ! -s "$result_root/preflight/target_contigs_missing_from_fasta.txt" && ! -s "$result_root/preflight/background_contigs_missing_from_fasta.txt" ]] || { printf 'BED/FASTA contig incompatibility.\n' >&2; exit 1; }

{
    printf 'metric\tvalue\n'
    printf 'target_regions\t%s\n' "$target_n"
    printf 'background_regions\t%s\n' "$background_n"
    printf 'target_unique_contigs\t%s\n' "$(wc -l < "$result_root/preflight/target_contigs.txt")"
    printf 'background_unique_contigs\t%s\n' "$(wc -l < "$result_root/preflight/background_contigs.txt")"
    printf 'fasta_contigs\t%s\n' "$(wc -l < "$result_root/preflight/fasta_contigs.txt")"
    printf 'target_ids_missing_from_background\t0\n'
    printf 'target_contigs_missing_from_fasta\t0\n'
    printf 'background_contigs_missing_from_fasta\t0\n'
    printf 'central_window_bp\t200\n'
    printf 'processors\t1\n'
} > "$result_root/preflight/preflight_metrics.tsv"

{
    printf 'file\tsha256\n'
    for input in "$target_bed" "$background_bed" "$fasta_gz" "$homer_root/bin/findMotifsGenome.pl" "$homer_root/bin/HomerSVGLogo.pm" "$homer_root/bin/mergePeaks.exe" "$homer_root/bin/homerTools.exe" "$homer_root/bin/homer2.exe"; do
        printf '%s\t%s\n' "$input" "$(sha256sum "$input" | awk '{print $1}')"
    done
} > "$result_root/preflight/input_hashes.tsv"

{
    printf 'HOMER_PACKAGE_VERSION=5.1\n'
    printf 'HOMER_PACKAGE_SOURCE=bioconda_homer-5.1-pl5262h9948957_0.tar.bz2\n'
    printf 'HOMER_PACKAGE_SHA256=58436010bfebe77abc867204b629e93622df09abe51acbe41fe92815468c0570\n'
    printf 'CYGWIN_REBUILD=single_worker_gcc_14.4.0\n'
    printf 'PATCH_1=Archived 15 Linux binaries so Cygwin resolves the compiled .exe files.\n'
    printf 'PATCH_2=HomerSVGLogo.pm line 163 calls plotXY directly for Perl 5.44 compatibility.\n'
    printf 'PATCH_3=findMotifsGenome.pl skips mkdir when the output directory already exists.\n'
    printf 'PERL_VERSION=%s\n' "$(perl -e 'printf "%vd", $^V')"
    printf 'CYGWIN_VERSION=%s\n' "$(uname -srvmo)"
} > "$result_root/preflight/software_versions.txt"

{
    printf 'status\tPASS\n'
    printf 'scope\tExploratory nominal-candidate motif sensitivity analysis; not FDR-significant DAR evidence.\n'
    printf 'target_definition\tPrimary 10-in-5 raw p-value <0.05 plus absolute ordinary LFC >=1 (1,562 peaks).\n'
    printf 'background_definition\tAll 59,186 primary 10-in-5 tested peaks.\n'
    printf 'window\tFixed central 200 bp.\n'
    printf 'reference\tGENCODE Release 50 GRCh38 primary assembly FASTA; published MD5 verified.\n'
} > "$result_root/preflight/PREFLIGHT_STATUS.txt"

printf 'PREFLIGHT=PASS\n'
printf 'COMMAND=findMotifsGenome.pl %s %s %s -size 200 -bg %s -p 1\n' "$target_bed" "$fasta" "$result_root/homer_output" "$background_bed"
findMotifsGenome.pl "$target_bed" "$fasta" "$result_root/homer_output" -size 200 -bg "$background_bed" -p 1

for expected in "$result_root/homer_output/knownResults.txt" "$result_root/homer_output/knownResults.html" "$result_root/homer_output/homerResults.html"; do
    [[ -s "$expected" ]] || { printf 'Missing or empty expected HOMER output: %s\n' "$expected" >&2; exit 1; }
done

{
    printf 'WEEK3_HOMER_PUBLIC_REPRODUCTION\tPASS\n'
    printf 'target_regions\t%s\n' "$target_n"
    printf 'background_regions\t%s\n' "$background_n"
    printf 'central_window_bp\t200\n'
    printf 'processors\t1\n'
    printf 'interpretation\tExploratory nominal-candidate sensitivity analysis; zero FDR-significant DARs remains the primary result.\n'
} > "$result_root/STATUS.txt"

printf 'END_UTC=%s\n' "$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
printf 'WEEK3_HOMER_PUBLIC_REPRODUCTION=PASS\n'
