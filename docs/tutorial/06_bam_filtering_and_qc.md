# 06. BAM filtering and alignment QC

## Pair-aware processing

Convert, filter, coordinate-sort, add or repair read groups, handle duplicates,
and index BAMs in declared stages. A common starting rule is to retain properly
paired primary alignments above a mapping-quality threshold while excluding
unmapped, mate-unmapped, secondary, supplementary, and vendor-failed records.
The exact flag mask and MAPQ threshold must be reviewed for the study.

The internship used proper-pair filtering, MAPQ >= 30, and exclusion mask 1804.
This is documented evidence for that project, not a universal prescription.

## Duplicate and blacklist decisions

Record whether duplicates are marked or removed and why. For blacklist removal
in paired-end data, remove the entire fragment if either mate overlaps a
blacklisted interval. Filtering individual alignment records can leave orphaned
mates and must not be described as pair-preserving.

## Useful checks

```bash
samtools quickcheck -v sample.filtered.bam
samtools flagstat sample.filtered.bam > logs/sample.flagstat.txt
samtools idxstats sample.filtered.bam > logs/sample.idxstats.txt
```

Also review duplicate metrics, insert-size distribution, mitochondrial fraction
when relevant, nucleosome-associated fragment-length periodicity, TSS enrichment
against a compatible annotation, unexpected contigs, read-group identity,
coordinate sorting, BAM and BAI agreement, and pair counts before and after each
filter. Library-complexity metrics can help distinguish deep sequencing from
excessive duplication.

## BAM gate

The final BAM is nonempty, coordinate sorted, indexed, structurally valid, and
contains the intended sample only. Pair preservation and filtering counts are
explicitly validated before peak calling. Compare QC metrics across samples and
document any study-specific acceptance criteria before excluding a library.
