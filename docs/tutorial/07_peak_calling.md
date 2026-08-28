# 07. Peak calling

## Call peaks per biological sample

For paired-end ATAC-seq, MACS3 can infer fragments directly from paired
alignments with `-f BAMPE`. Confirm whether narrow or broad peaks fit the
biological question and keep the choice consistent across comparable samples.

```bash
macs3 callpeak \
  -t sample.final.bam \
  -f BAMPE \
  -g GENOME_SIZE \
  -n sample \
  --outdir results/peaks/sample \
  -q 0.05
```

`GENOME_SIZE`, the q-value, optional control inputs, duplicate handling, and
additional options require study-specific review. The example is not suitable
for every assay or organism.

## Peak gate

- the caller reports successful completion;
- output BED-like files have valid, sorted coordinates on expected contigs;
- peak counts, widths, and signal distributions are compared across samples;
- fraction of reads in peaks (FRiP) is calculated consistently after peaks are
  available and interpreted together with depth, TSS enrichment, and complexity;
- low-quality or unusual samples are investigated alongside alignment QC;
- peak files and the exact BAM inputs are recorded with checksums;
- biological replicates remain separate at this stage.
