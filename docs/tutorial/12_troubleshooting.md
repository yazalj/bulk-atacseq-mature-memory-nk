# 12. Troubleshooting and stop conditions

| Observation | Investigate | Do not do automatically |
|---|---|---|
| Missing or unequal R1/R2 records | transfer integrity, truncation, lane composition, checksums | continue with orphaned paired-end inputs |
| Unexpectedly high adapter content | library protocol, read length, adapter identity | guess adapter sequences |
| Low or uneven mapping | reference build, contamination, read quality, aligner settings | exclude a sample without documented review |
| BAM validation failure | truncation, sort order, header, index freshness | call peaks from the failed BAM |
| One mate remains after filtering | record-level rather than fragment-level filtering | describe the output as pair-preserving |
| Large peak-count differences | library depth, duplicate rate, TSS/fragment QC, peak parameters | merge away the evidence before review |
| Count columns and metadata differ | basenames, explicit ordering, renamed samples | reorder by visual inspection alone |
| Rank-deficient model | confounding, empty factor combinations, duplicated covariates | remove terms until software runs without scientific justification |
| Extreme size factor or PCA outlier | depth, sample identity, QC history, batch, contamination | delete the sample because of one plot |
| Zero FDR-significant regions | power, heterogeneity, model, diagnostics, effect sizes | promote raw-p-value candidates to discoveries |
| Enrichment changes with background | target universe and gene mapping | use the whole genome by convenience |

## Copy-and-adapt diagnostic commands

Use only the checks relevant to the failed stage. Replace the paths before
running them; these commands inspect inputs and do not repair or replace them.

```bash
# Compressed FASTQ integrity
gzip -t /path/to/sample_R1.fastq.gz
gzip -t /path/to/sample_R2.fastq.gz

# Read counts: each FASTQ record contains four lines
zcat /path/to/sample_R1.fastq.gz | awk 'END { print NR / 4 }'
zcat /path/to/sample_R2.fastq.gz | awk 'END { print NR / 4 }'

# BAM structure, header, counts, and index agreement
samtools quickcheck -v /path/to/sample.final.bam
samtools view -H /path/to/sample.final.bam
samtools flagstat /path/to/sample.final.bam
samtools idxstats /path/to/sample.final.bam

# Confirm the count-matrix header and metadata order before opening R
head -n 1 results/counts/consensus_counts.clean.tsv
cut -f 2 config/local/samples.tsv
```

If one of these checks fails, preserve the failing file and log, determine
whether the problem arose during download, transfer, processing, or indexing,
and create a new attempt rather than overwriting the evidence.

## Stop and resolve before continuing when

- a sample identity, mate, organism, reference build, annotation, blacklist, or
  contrast is uncertain;
- an input is truncated, writable when it should be immutable, or fails its
  checksum;
- a stage would overwrite an existing result;
- coordinate conventions or chromosome names disagree;
- BAM structure, sample count, pair preservation, matrix order, or model rank
  fails validation;
- compute limits or tool behavior are unclear;
- the next interpretation would require changing the prespecified endpoint.

Record failures, their causes, and the chosen resolution. A stopped workflow
with preserved evidence is more reproducible than a completed workflow built on
an unresolved input or design problem.
