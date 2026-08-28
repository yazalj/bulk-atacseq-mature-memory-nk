# 01. Requirements and project setup

## Recommended environment

Most command-line ATAC-seq tools are easiest to run on Linux or a managed HPC
environment. Windows users can use a validated Linux environment such as WSL
when permitted by their institution. The R differential-analysis stage can be
run on any supported platform if file formats and sample order remain stable.

Core tools used in the worked example were FastQC, MultiQC, Cutadapt, Bowtie2,
SAMtools, Picard, bedtools, MACS3, featureCounts, R, DESeq2, and apeglm. See the
[recorded case-study versions](../../environment/software_versions.tsv).
Version compatibility should be checked for a new run; matching a version
number does not replace validation.

## Suggested working layout

```text
project/
├── config/       metadata, input manifest, and parameters
├── workflow/     scripts kept under version control
├── raw/          immutable reads, normally outside Git
├── reference/    FASTA, indexes, annotation, and blacklist, outside Git
├── results/      stage-specific, non-overwriting outputs
├── logs/         commands, versions, warnings, and validation summaries
└── docs/         decisions, methods, limitations, and report material
```

## Setup gate

Before analysis, confirm that:

- each required executable resolves to the intended installation;
- versions are recorded;
- compute limits are known and explicitly set;
- raw and reference inputs are read-only;
- outputs cannot overwrite inputs or previous attempts;
- available disk space is suitable for FASTQ, BAM, and temporary files;
- secrets, controlled data, and absolute private paths will not enter Git.

Copy the three templates from [`config/templates/`](../../config/templates/)
into a new project-specific configuration area, then replace every placeholder.
