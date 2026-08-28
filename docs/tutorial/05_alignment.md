# 05. Alignment

## Prepare the reference

Use an authoritative FASTA and an aligner index built from exactly that FASTA.
Record the organism, assembly, source, release, checksum, included contigs, and
chromosome naming convention. Annotation and blacklist resources used later
must be compatible with the same assembly.

## Paired-end Bowtie2 pattern

```bash
bowtie2 \
  --threads 2 \
  -x /path/to/index_prefix \
  -1 results/trimmed/sample_R1.trimmed.fastq.gz \
  -2 results/trimmed/sample_R2.trimmed.fastq.gz \
  -S results/alignment/sample.sam \
  2> logs/sample.bowtie2.log
```

Bowtie2 alignment presets and fragment-length assumptions should match the
library. Do not treat an example command as a validated protocol for a new
dataset.

## Alignment gate

- exit status is zero and the log is preserved;
- the input read-pair count agrees with the alignment summary;
- headers identify the expected reference;
- overall and concordant mapping rates are reviewed across all samples;
- no sample is accepted or excluded solely by comparison with an undocumented
  threshold.
