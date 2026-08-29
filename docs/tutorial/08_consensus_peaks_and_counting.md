# 08. Consensus peaks and fragment counting

[← Previous: Peak calling](07_peak_calling.md) · [Tutorial contents](README.md) · [Next: Differential accessibility →](09_differential_accessibility.md)

## Define the feature universe

A differential analysis needs the same genomic features for every sample.
Declare how replicate peak sets become a consensus universe. Options include a
simple union, a reproducibility-supported set, or a study-specific consensus
rule. The choice affects both sensitivity and the multiple-testing burden.

The internship merged peaks within each condition, combined the condition
unions, and merged overlapping or book-ended intervals with `bedtools merge
-d 0`. Singleton-supported regions were retained. This created 112,759 regions.

The worked implementation is
[`workflow/week2_consensus_counts/`](../../workflow/week2_consensus_counts/).
Its safety checks are intentionally locked to the internship's ten samples, so
a new study must revise the manifest and expected dimensions deliberately.

## Copy-and-adapt consensus commands

The following manual example creates a union within each of two conditions and
then merges both condition unions. Replace the peak arrays with every biological
sample in the corresponding condition. Do not use a wildcard that might include
an unrelated or failed peak file.

```bash
set -euo pipefail

CONDITION_A_PEAKS=(
  "results/peaks/sample_01/sample_01_peaks.narrowPeak"
  "results/peaks/sample_03/sample_03_peaks.narrowPeak"
)

CONDITION_B_PEAKS=(
  "results/peaks/sample_02/sample_02_peaks.narrowPeak"
  "results/peaks/sample_04/sample_04_peaks.narrowPeak"
)

mkdir -p results/consensus results/counts logs/counts
: > results/consensus/condition_A.all_peaks.bed
: > results/consensus/condition_B.all_peaks.bed

for peak_file in "${CONDITION_A_PEAKS[@]}"; do
  test -s "${peak_file}"
  cut -f 1-3 "${peak_file}" >> results/consensus/condition_A.all_peaks.bed
done

for peak_file in "${CONDITION_B_PEAKS[@]}"; do
  test -s "${peak_file}"
  cut -f 1-3 "${peak_file}" >> results/consensus/condition_B.all_peaks.bed
done

sort -k1,1 -k2,2n -k3,3n \
  results/consensus/condition_A.all_peaks.bed \
  > results/consensus/condition_A.sorted.bed
bedtools merge -d 0 \
  -i results/consensus/condition_A.sorted.bed \
  > results/consensus/condition_A.union.bed

sort -k1,1 -k2,2n -k3,3n \
  results/consensus/condition_B.all_peaks.bed \
  > results/consensus/condition_B.sorted.bed
bedtools merge -d 0 \
  -i results/consensus/condition_B.sorted.bed \
  > results/consensus/condition_B.union.bed

cat results/consensus/condition_A.union.bed \
    results/consensus/condition_B.union.bed \
  > results/consensus/all_conditions.bed

sort -k1,1 -k2,2n -k3,3n \
  results/consensus/all_conditions.bed \
  > results/consensus/all_conditions.sorted.bed
bedtools merge -d 0 \
  -i results/consensus/all_conditions.sorted.bed \
  > results/consensus/consensus.bed

awk 'BEGIN { OFS="\t"; print "GeneID", "Chr", "Start", "End", "Strand" }
     { printf "peak_%06d\t%s\t%d\t%d\t.\n", NR, $1, $2 + 1, $3 }' \
  results/consensus/consensus.bed \
  > results/consensus/consensus.saf

wc -l results/consensus/*.union.bed results/consensus/consensus.bed
```

This is a union strategy, matching the internship's teaching approach. If the
study requires replicate support or another consensus definition, replace this
block and document the alternative before counting.

## Count paired fragments

Convert the BED universe to the annotation format required by the counter while
preserving coordinate semantics. The internship used one-based inclusive SAF
coordinates and featureCounts:

```text
-F SAF -p --countReadPairs -B -C -s 0 -T 2
```

Copy the BAMs from `config/local/samples.tsv` into the array in exactly the
same order as the metadata rows:

```bash
set -euo pipefail

BAMS=(
  "results/bam/sample_01.final.bam"
  "results/bam/sample_02.final.bam"
  "results/bam/sample_03.final.bam"
  "results/bam/sample_04.final.bam"
)

for bam in "${BAMS[@]}"; do
  test -s "${bam}"
  test -s "${bam}.bai"
  samtools quickcheck -v "${bam}"
done

featureCounts \
  -F SAF \
  -a results/consensus/consensus.saf \
  -o results/counts/consensus_featurecounts.tsv \
  -p \
  --countReadPairs \
  -B \
  -C \
  -s 0 \
  -T 2 \
  "${BAMS[@]}" \
  2>&1 | tee logs/counts/featurecounts.log
```

The example is unstranded and requires both mates to be aligned to the same
chromosome. Review those choices for the new libraries.

## Create a clean count matrix without changing the counts

Paste this block into R after filling `config/local/samples.tsv`. It verifies
the featureCounts BAM order against the `final_bam` column before replacing BAM
column labels with stable sample IDs:

```r
metadata <- read.delim(
  "config/local/samples.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

featurecounts <- read.delim(
  "results/counts/consensus_featurecounts.tsv",
  comment.char = "#",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (ncol(featurecounts) != 6L + nrow(metadata)) {
  stop("featureCounts sample-column number does not match metadata rows.")
}

observed_bams <- basename(colnames(featurecounts)[7:ncol(featurecounts)])
expected_bams <- basename(metadata$final_bam)
if (!identical(observed_bams, expected_bams)) {
  stop("featureCounts BAM order does not exactly match metadata final_bam order.")
}

count_matrix <- featurecounts[, 7:ncol(featurecounts), drop = FALSE]
colnames(count_matrix) <- metadata$sample_id
rownames(count_matrix) <- featurecounts$Geneid

if (anyNA(count_matrix) || any(as.matrix(count_matrix) < 0) ||
    any(as.matrix(count_matrix) %% 1 != 0)) {
  stop("Counts must be nonnegative integers without missing values.")
}

clean_output <- data.frame(
  peak_id = rownames(count_matrix),
  count_matrix,
  check.names = FALSE
)

write.table(
  clean_output,
  "results/counts/consensus_counts.clean.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  featurecounts[, c("Geneid", "Chr", "Start", "End", "Strand", "Length")],
  "results/counts/consensus_peak_annotation.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
```

## Counting gate

- consensus intervals are valid, sorted, and non-overlapping;
- coordinate conversion preserves interval lengths exactly;
- BAM order is explicit and matches the planned metadata order;
- counts are nonnegative integers and represent fragments, not individual mates;
- column sums agree with the counter's assignment summary;
- the clean matrix contains no silent filtering, normalization, or reordering.

---

[← Previous: Peak calling](07_peak_calling.md) · [Tutorial contents](README.md) · [Next: Differential accessibility →](09_differential_accessibility.md)
