# Resource planning before real data

[← Previous: Overview](00_overview.md) · [Tutorial contents](README.md) · [Next: Requirements and setup →](01_requirements_and_project_setup.md)

Bulk ATAC-seq creates large intermediate files. Check storage, memory, runtime,
and institutional computing rules before downloading reads. The ranges below
are conservative planning aids, not performance guarantees or benchmark
claims.

## Planning assumptions

The table assumes a human paired-end library with roughly 30–50 million read
pairs, conventional compressed FASTQ input, and one intensive job at a time.
The tutorial commands default to two total CPU threads. Reference size, read
length, sequencing depth, filesystem speed, compression, and tool versions can
move a run far outside these ranges. Pilot one representative sample before
estimating the complete cohort.

| Stage | Rough elapsed time | Rough peak memory | Additional working disk | Main scaling factor |
|---|---:|---:|---:|---|
| Environment installation | 15–90 minutes once | 2–8 GB | 8–20 GB | package solving and download speed |
| Reference download and Bowtie2 index | 1–6 hours once | 8–32 GB | 15–40 GB | genome size and index construction |
| Raw FastQC and focused MultiQC | 10–45 minutes per sample | 2–4 GB | usually below 1 GB | read count and compression speed |
| Paired trimming and post-trim QC | 20–120 minutes per sample | 2–8 GB | about one additional FASTQ copy | read count and compression speed |
| Paired alignment | 1–8 hours per sample | 8–32 GB | roughly 1–3 times the compressed FASTQ size | read count, genome, and disk speed |
| BAM filtering, sorting, duplicate handling, and QC | 1–6 hours per sample | 8–32 GB | temporarily 2–4 times one BAM | BAM size and sorting temporary space |
| Peak calling | 10–90 minutes per sample | 4–16 GB | commonly below 2 GB | usable fragment count |
| Consensus construction and fragment counting | 1–8 hours for a small cohort | 8–32 GB | commonly 1–10 GB | sample and consensus-region counts |
| Differential accessibility | 10 minutes–3 hours | 8–32 GB | commonly below 5 GB | tested regions, samples, and model complexity |
| Annotation, enrichment, motifs, and locus review | 15 minutes–many hours | 4–32 GB | roughly 1–20 GB | target/background size and motif search scope |

These numbers describe ordinary teaching-scale work. Deep cohorts, large
genomes, many conditions, or broad motif searches can require much more.

## Storage rule of thumb

Record the total size of the compressed FASTQ files, then plan for at least
three to five times that amount as working space when FASTQ, BAM, sorted BAM,
temporary sorting files, indexes, and retained checkpoints coexist. More is
needed when keeping multiple attempts. Never use the final few gigabytes of a
filesystem for BAM sorting or index creation.

From the intended project location, inspect the available filesystem:

```bash
df -h .
```

After inputs are present, inspect the project footprint without changing it:

```bash
du -sh .
du -sh data results logs 2>/dev/null || true
```

## Before each stage

1. Identify the exact inputs and estimate their combined size.
2. Confirm free disk space and the location used for temporary files.
3. Confirm the scheduler or workstation memory and thread limits.
4. Run one sample or a small subset first when the stage supports it.
5. Inspect its output size, elapsed time, log, and validation evidence.
6. Recalculate the cohort estimate before continuing.

Do not launch several CPU- or disk-intensive samples simultaneously merely to
save wall-clock time. On managed infrastructure, follow the administrator's
limits even when the machine appears idle.

## Fast technical confidence check

Before downloading real data, run the repository's dependency-free
[synthetic smoke test](../../demo/README.md). It checks ordered metadata and
count-matrix validation in less than a second, but it does not estimate the
runtime of sequencing tools.

---

[← Previous: Overview](00_overview.md) · [Tutorial contents](README.md) · [Next: Requirements and setup →](01_requirements_and_project_setup.md)
