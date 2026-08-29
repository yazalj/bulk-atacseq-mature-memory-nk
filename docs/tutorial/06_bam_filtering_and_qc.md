# 06. BAM filtering and alignment QC

## Pair-aware processing

Convert, filter, coordinate-sort, add or repair read groups, handle duplicates,
and index BAMs in declared stages. A common starting rule is to retain properly
paired primary alignments above a mapping-quality threshold while excluding
unmapped, mate-unmapped, secondary, supplementary, and vendor-failed records.
The exact flag mask and MAPQ threshold must be reviewed for the study.

The internship used proper-pair filtering, MAPQ >= 30, and exclusion mask 1804.
This is documented evidence for that project, not a universal prescription.

## Copy-and-adapt filtering and sorting commands

The block below uses the internship rule as a visible example. Change
`MIN_MAPQ`, `REQUIRE_FLAGS`, and `EXCLUDE_FLAGS` only after recording the rule
for the new study. It records every QNAME for which any alignment fails the
rule, then removes the entire QNAME. This prevents record-level MAPQ or flag
filtering from leaving one retained mate. Separate commands are used instead of
one long pipe so every intermediate can be checked.

```bash
set -euo pipefail

SAMPLE_ID="sample_01"
INPUT_SAM="results/alignment/${SAMPLE_ID}.sam"
MIN_MAPQ=30
REQUIRE_FLAGS=2
EXCLUDE_FLAGS=1804

mkdir -p results/bam logs/bam

samtools view \
  -@ 1 \
  -b \
  -q "${MIN_MAPQ}" \
  -f "${REQUIRE_FLAGS}" \
  -F "${EXCLUDE_FLAGS}" \
  -o "results/bam/${SAMPLE_ID}.passing_records.bam" \
  -U "results/bam/${SAMPLE_ID}.failing_records.bam" \
  "${INPUT_SAM}"

samtools view "results/bam/${SAMPLE_ID}.failing_records.bam" \
  | cut -f 1 \
  | sort -u \
  > "results/bam/${SAMPLE_ID}.failing_qnames.txt"

if [[ -s "results/bam/${SAMPLE_ID}.failing_qnames.txt" ]]; then
  samtools view \
    -b \
    -N "^results/bam/${SAMPLE_ID}.failing_qnames.txt" \
    -o "results/bam/${SAMPLE_ID}.filtered.unsorted.bam" \
    "${INPUT_SAM}"
else
  samtools view \
    -b \
    -o "results/bam/${SAMPLE_ID}.filtered.unsorted.bam" \
    "${INPUT_SAM}"
fi

samtools sort \
  -@ 1 \
  -o "results/bam/${SAMPLE_ID}.filtered.sorted.bam" \
  "results/bam/${SAMPLE_ID}.filtered.unsorted.bam"

samtools quickcheck -v "results/bam/${SAMPLE_ID}.filtered.sorted.bam"
samtools flagstat "results/bam/${SAMPLE_ID}.filtered.sorted.bam" \
  > "logs/bam/${SAMPLE_ID}.filtered.flagstat.txt"

samtools view "results/bam/${SAMPLE_ID}.filtered.sorted.bam" \
  | cut -f 1 \
  | sort \
  | uniq -c \
  | awk '$1 != 2 { print }' \
  > "logs/bam/${SAMPLE_ID}.unexpected_qname_multiplicity.txt"

test ! -s "logs/bam/${SAMPLE_ID}.unexpected_qname_multiplicity.txt"
echo "Pair-aware MAPQ and flag filtering: PASS"
```

## Add or confirm read groups

First inspect the header:

```bash
samtools view -H results/bam/sample_01.filtered.sorted.bam | grep '^@RG' || true
```

If a correct sample-level read group is already present, do not replace it. If
it is missing, adapt and run:

```bash
SAMPLE_ID="sample_01"
READ_GROUP_PLATFORM_UNIT="CHANGE_ME_LIBRARY_OR_LANE"

if [[ "${READ_GROUP_PLATFORM_UNIT}" == "CHANGE_ME_LIBRARY_OR_LANE" ]]; then
  echo "Replace the read-group platform unit before continuing." >&2
  exit 1
fi

picard AddOrReplaceReadGroups \
  I="results/bam/${SAMPLE_ID}.filtered.sorted.bam" \
  O="results/bam/${SAMPLE_ID}.filtered.sorted.rg.bam" \
  RGID="${SAMPLE_ID}" \
  RGLB="${SAMPLE_ID}" \
  RGPL="ILLUMINA" \
  RGPU="${READ_GROUP_PLATFORM_UNIT}" \
  RGSM="${SAMPLE_ID}" \
  VALIDATION_STRINGENCY=STRICT
```

## Duplicate handling

The following example physically removes duplicates. If the prespecified policy
is to mark and retain them, set `REMOVE_DUPLICATES=false` and account for that
choice during peak calling.

```bash
SAMPLE_ID="sample_01"
INPUT_RG_BAM="results/bam/${SAMPLE_ID}.filtered.sorted.rg.bam"

picard MarkDuplicates \
  I="${INPUT_RG_BAM}" \
  O="results/bam/${SAMPLE_ID}.deduplicated.bam" \
  M="logs/bam/${SAMPLE_ID}.duplicate_metrics.txt" \
  REMOVE_DUPLICATES=true \
  CREATE_INDEX=true \
  VALIDATION_STRINGENCY=STRICT
```

## Duplicate and blacklist decisions

Record whether duplicates are marked or removed and why. For blacklist removal
in paired-end data, remove the entire fragment if either mate overlaps a
blacklisted interval. Filtering individual alignment records can leave orphaned
mates and must not be described as pair-preserving.

If a compatible blacklist is part of the declared method, this block first
collects every overlapping read name and then removes **all** records with any
of those names. This preserves the paired-fragment decision even when only one
mate overlaps. It requires a recent SAMtools version that supports negated
`-N` read-name lists.

```bash
set -euo pipefail

SAMPLE_ID="sample_01"
INPUT_BAM="results/bam/${SAMPLE_ID}.deduplicated.bam"
BLACKLIST_BED="/path/to/build-compatible.blacklist.bed"
OVERLAP_BAM="results/bam/${SAMPLE_ID}.blacklist_overlaps.bam"
BLACKLIST_QNAMES="results/bam/${SAMPLE_ID}.blacklist_qnames.txt"
FINAL_BAM="results/bam/${SAMPLE_ID}.final.bam"

test -s "${INPUT_BAM}"
test -s "${BLACKLIST_BED}"

bedtools intersect \
  -abam "${INPUT_BAM}" \
  -b "${BLACKLIST_BED}" \
  -u \
  > "${OVERLAP_BAM}"

samtools view "${OVERLAP_BAM}" | cut -f 1 | sort -u > "${BLACKLIST_QNAMES}"

if [[ -s "${BLACKLIST_QNAMES}" ]]; then
  samtools view \
    -b \
    -N "^${BLACKLIST_QNAMES}" \
    -o "${FINAL_BAM}" \
    "${INPUT_BAM}"
else
  cp "${INPUT_BAM}" "${FINAL_BAM}"
fi

samtools index "${FINAL_BAM}"
samtools quickcheck -v "${FINAL_BAM}"

bedtools intersect \
  -abam "${FINAL_BAM}" \
  -b "${BLACKLIST_BED}" \
  -u \
  > "results/bam/${SAMPLE_ID}.post_blacklist_overlap_check.bam"

test "$(samtools view -c "results/bam/${SAMPLE_ID}.post_blacklist_overlap_check.bam")" -eq 0
echo "Pair-preserving blacklist validation: PASS"
```

If no suitable blacklist exists for the organism and assembly, do not use an
unrelated one. Record that the stage was intentionally omitted.

## Useful checks

```bash
samtools quickcheck -v sample.filtered.bam
samtools flagstat sample.filtered.bam > logs/sample.flagstat.txt
samtools idxstats sample.filtered.bam > logs/sample.idxstats.txt
```

Collect insert-size metrics from the final paired BAM:

```bash
SAMPLE_ID="sample_01"

picard CollectInsertSizeMetrics \
  I="results/bam/${SAMPLE_ID}.final.bam" \
  O="logs/bam/${SAMPLE_ID}.insert_size_metrics.txt" \
  H="results/bam/${SAMPLE_ID}.insert_size_histogram.pdf" \
  VALIDATION_STRINGENCY=STRICT
```

Also review duplicate metrics, insert-size distribution, mitochondrial fraction
when relevant, nucleosome-associated fragment-length periodicity, TSS enrichment
against a compatible annotation, unexpected contigs, read-group identity,
coordinate sorting, BAM and BAI agreement, and pair counts before and after each
filter. Library-complexity metrics can help distinguish deep sequencing from
excessive duplication.

## BAM gate

The final BAM is nonempty, coordinate sorted, indexed, structurally valid, and
contains the intended sample only. Pair preservation and filtering counts are
explicitly validated before peak calling. Compare QC metrics across samples and
document any study-specific acceptance criteria before excluding a library.
