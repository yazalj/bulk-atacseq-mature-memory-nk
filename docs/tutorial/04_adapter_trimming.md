# 04. Adapter and quality trimming

## Decide from evidence

Trim only when adapters or unusable terminal bases are supported by the library
preparation and raw-read QC. Confirm the actual adapter sequences; do not copy
them from an unrelated protocol.

## Paired-end command pattern

```bash
cutadapt \
  -a ADAPTER_READ1 \
  -A ADAPTER_READ2 \
  -o results/trimmed/sample_R1.trimmed.fastq.gz \
  -p results/trimmed/sample_R2.trimmed.fastq.gz \
  sample_R1.fastq.gz sample_R2.fastq.gz
```

Add minimum-length or quality options only when justified and record them.
Always process the two mates together so pair correspondence is preserved.

## Validation gate

- Cutadapt exits successfully and reports the expected input pairs;
- output mates are nonempty and have equal numbers of records;
- trimming rates are plausible and not dominated by unexpectedly short reads;
- FastQC is rerun on the trimmed mates;
- raw reads remain unchanged.
