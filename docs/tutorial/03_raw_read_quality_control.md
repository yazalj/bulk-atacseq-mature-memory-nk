# 03. Raw-read quality control

[← Previous: Samples and metadata](02_samples_and_metadata.md) · [Tutorial contents](README.md) · [Next: Adapter trimming →](04_adapter_trimming.md)

## Purpose

Inspect read quality, length, base composition, duplication, overrepresented
sequences, and possible adapter content before changing the data.

## Copy-and-adapt commands for one sample

Replace the three values below, then paste the block from the project root.
Repeat the FastQC command for every sample before running MultiQC.

```bash
set -euo pipefail

SAMPLE_ID="sample_01"
R1="/path/to/sample_01_R1.fastq.gz"
R2="/path/to/sample_01_R2.fastq.gz"

test -s "${R1}"
test -s "${R2}"
mkdir -p results/qc/raw logs/qc

fastqc \
  --threads 2 \
  --outdir results/qc/raw \
  "${R1}" "${R2}" \
  2>&1 | tee "logs/qc/${SAMPLE_ID}.raw_fastqc.log"

# Run this summary after FastQC has completed for every raw mate.
multiqc \
  --outdir results/qc/raw_multiqc \
  results/qc/raw \
  2>&1 | tee logs/qc/raw_multiqc.log
```

Run FastQC for every mate. Restrict MultiQC to the intended stage so raw and
trimmed reports are not silently mixed.

Confirm that each sample produced one R1 and one R2 HTML report and one R1 and
one R2 ZIP archive before deciding whether trimming is needed.

## Interpretation

ATAC-seq libraries may show non-random sequence composition and duplication
that reflect the assay rather than a technical failure. Adapter evidence,
quality decay, abnormal read lengths, severe mate imbalance, or a sample that
differs sharply from its peers should be investigated in context.

## QC gate

Record the number of files, read lengths, read counts, adapter evidence,
important warnings, and the decision to trim or retain reads. Do not discard a
sample solely because a FastQC module is colored red.

---

[← Previous: Samples and metadata](02_samples_and_metadata.md) · [Tutorial contents](README.md) · [Next: Adapter trimming →](04_adapter_trimming.md)
