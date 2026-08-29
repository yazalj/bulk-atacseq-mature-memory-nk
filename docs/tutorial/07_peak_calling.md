# 07. Peak calling

## Call peaks per biological sample

For paired-end ATAC-seq, MACS3 can infer fragments directly from paired
alignments with `-f BAMPE`. Confirm whether narrow or broad peaks fit the
biological question and keep the choice consistent across comparable samples.

```bash
set -euo pipefail

SAMPLE_ID="sample_01"
FINAL_BAM="results/bam/${SAMPLE_ID}.final.bam"
GENOME_SIZE="CHANGE_ME"
QVALUE=0.05

if [[ "${GENOME_SIZE}" == "CHANGE_ME" ]]; then
  echo "Replace GENOME_SIZE for the selected organism and assembly." >&2
  exit 1
fi

test -s "${FINAL_BAM}"
mkdir -p "results/peaks/${SAMPLE_ID}" logs/peaks

macs3 callpeak \
  -t "${FINAL_BAM}" \
  -f BAMPE \
  -g "${GENOME_SIZE}" \
  -n "${SAMPLE_ID}" \
  --outdir "results/peaks/${SAMPLE_ID}" \
  -q "${QVALUE}" \
  --keep-dup all \
  2>&1 | tee "logs/peaks/${SAMPLE_ID}.macs3.log"

test -s "results/peaks/${SAMPLE_ID}/${SAMPLE_ID}_peaks.narrowPeak"
```

`GENOME_SIZE`, the q-value, optional control inputs, duplicate handling, and
additional options require study-specific review. The example is not suitable
for every assay or organism. It uses `--keep-dup all` because Chapter 6's
example physically removes duplicates before peak calling. If duplicates were
retained instead, review the MACS3 duplicate policy rather than copying this
option automatically.

## Peak gate

- the caller reports successful completion;
- output BED-like files have valid, sorted coordinates on expected contigs;
- peak counts, widths, and signal distributions are compared across samples;
- fraction of reads in peaks (FRiP) is calculated consistently after peaks are
  available and interpreted together with depth, TSS enrichment, and complexity;
- low-quality or unusual samples are investigated alongside alignment QC;
- peak files and the exact BAM inputs are recorded with checksums;
- biological replicates remain separate at this stage.
