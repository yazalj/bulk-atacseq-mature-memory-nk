# 05. Alignment

## Prepare the reference

Use an authoritative FASTA and an aligner index built from exactly that FASTA.
Record the organism, assembly, source, release, checksum, included contigs, and
chromosome naming convention. Annotation and blacklist resources used later
must be compatible with the same assembly.

Build the Bowtie2 index once for a newly obtained reference. Do not rebuild it
for every sample:

```bash
set -euo pipefail

REFERENCE_FASTA="/path/to/reference.fa"
INDEX_PREFIX="reference/bowtie2/genome"

test -s "${REFERENCE_FASTA}"
mkdir -p reference/bowtie2 logs/reference

bowtie2-build \
  --threads 2 \
  "${REFERENCE_FASTA}" "${INDEX_PREFIX}" \
  2>&1 | tee logs/reference/bowtie2_build.log
```

## Copy-and-adapt alignment commands for one sample

```bash
set -euo pipefail

SAMPLE_ID="sample_01"
INDEX_PREFIX="reference/bowtie2/genome"
R1="results/trimmed/${SAMPLE_ID}_R1.trimmed.fastq.gz"
R2="results/trimmed/${SAMPLE_ID}_R2.trimmed.fastq.gz"

test -s "${R1}"
test -s "${R2}"
mkdir -p results/alignment logs/alignment

bowtie2 \
  --threads 2 \
  -x "${INDEX_PREFIX}" \
  -1 "${R1}" \
  -2 "${R2}" \
  -S "results/alignment/${SAMPLE_ID}.sam" \
  2> "logs/alignment/${SAMPLE_ID}.bowtie2.log"

test -s "results/alignment/${SAMPLE_ID}.sam"
```

Bowtie2 alignment presets and fragment-length assumptions should match the
library. Do not treat an example command as a validated protocol for a new
dataset. If reads were intentionally left untrimmed, change `R1` and `R2` to
the validated raw FASTQ paths.

## Alignment gate

- exit status is zero and the log is preserved;
- the input read-pair count agrees with the alignment summary;
- headers identify the expected reference;
- overall and concordant mapping rates are reviewed across all samples;
- no sample is accepted or excluded solely by comparison with an undocumented
  threshold.
