# 04. Adapter and quality trimming

## Decide from evidence

Trim only when adapters or unusable terminal bases are supported by the library
preparation and raw-read QC. Confirm the actual adapter sequences; do not copy
them from an unrelated protocol.

## Copy-and-adapt paired-end commands

Use the adapter sequences documented for the actual library preparation. The
two example values deliberately refuse to run until they are changed.

```bash
set -euo pipefail

SAMPLE_ID="sample_01"
R1="/path/to/sample_01_R1.fastq.gz"
R2="/path/to/sample_01_R2.fastq.gz"
ADAPTER_R1="CHANGE_ME"
ADAPTER_R2="CHANGE_ME"

if [[ "${ADAPTER_R1}" == "CHANGE_ME" || "${ADAPTER_R2}" == "CHANGE_ME" ]]; then
  echo "Replace both adapter sequences before trimming." >&2
  exit 1
fi

mkdir -p results/trimmed results/qc/trimmed logs/trimming logs/qc

cutadapt \
  --cores 2 \
  -a "${ADAPTER_R1}" \
  -A "${ADAPTER_R2}" \
  -o "results/trimmed/${SAMPLE_ID}_R1.trimmed.fastq.gz" \
  -p "results/trimmed/${SAMPLE_ID}_R2.trimmed.fastq.gz" \
  "${R1}" "${R2}" \
  2>&1 | tee "logs/trimming/${SAMPLE_ID}.cutadapt.log"

fastqc \
  --threads 2 \
  --outdir results/qc/trimmed \
  "results/trimmed/${SAMPLE_ID}_R1.trimmed.fastq.gz" \
  "results/trimmed/${SAMPLE_ID}_R2.trimmed.fastq.gz" \
  2>&1 | tee "logs/qc/${SAMPLE_ID}.trimmed_fastqc.log"
```

Add minimum-length or quality options only when justified and record them.
Always process the two mates together so pair correspondence is preserved.
If the raw-read review does not justify trimming, record that decision and use
the raw mates as the alignment inputs instead of running this chapter.

## Validation gate

- Cutadapt exits successfully and reports the expected input pairs;
- output mates are nonempty and have equal numbers of records;
- trimming rates are plausible and not dominated by unexpectedly short reads;
- FastQC is rerun on the trimmed mates;
- raw reads remain unchanged.
