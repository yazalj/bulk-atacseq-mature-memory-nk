# 03. Raw-read quality control

## Purpose

Inspect read quality, length, base composition, duplication, overrepresented
sequences, and possible adapter content before changing the data.

## Example pattern

```bash
fastqc --threads 2 --outdir results/qc/raw sample_R1.fastq.gz sample_R2.fastq.gz
multiqc --outdir results/qc/raw_multiqc results/qc/raw
```

Run FastQC for every mate. Restrict MultiQC to the intended stage so raw and
trimmed reports are not silently mixed.

## Interpretation

ATAC-seq libraries may show non-random sequence composition and duplication
that reflect the assay rather than a technical failure. Adapter evidence,
quality decay, abnormal read lengths, severe mate imbalance, or a sample that
differs sharply from its peers should be investigated in context.

## QC gate

Record the number of files, read lengths, read counts, adapter evidence,
important warnings, and the decision to trim or retain reads. Do not discard a
sample solely because a FastQC module is colored red.
