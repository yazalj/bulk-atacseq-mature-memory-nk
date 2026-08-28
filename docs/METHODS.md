# Methods

## Week 1 training workflow

SRR7650763 is an untreated Immature NK paired-end sample analyzed on a GRCh38
chromosome-22 teaching subset. It was not included in the Mature versus Memory
comparison. FastQC assessed raw and trimmed reads. Cutadapt removed adapter
sequence, Bowtie2 aligned reads, and SAMtools retained properly paired records
with mapping quality at least 30 while excluding flag mask 1804. After
coordinate sorting and read-group repair, Picard MarkDuplicates removed
duplicate fragments. Blacklist removal was pair-preserving: if either mate
overlapped a blacklist interval, both alignments were removed.

MACS3 3.0.2 called peaks from the validated post-blacklist paired BAM in BAMPE
narrow-peak mode at q = 0.05. This chromosome-22 exercise is technical training
and does not provide a Mature-versus-Memory result.

## Consensus regions and fragment counting

The final comparison used four Mature NK and six Memory NK public runs. Peaks
were merged within each condition using `bedtools merge -d 0`. The condition
unions were then combined and merged with the same rule, retaining singleton-
supported regions. The resulting universe contained 112,759 regions.

The BED universe was converted to SAF with one-based inclusive starts.
featureCounts 2.1.1 counted paired fragments with:

```text
-F SAF -p --countReadPairs -B -C -s 0 -T 2
```

The cleaning step separated the raw featureCounts table into a count matrix and
peak annotation without filtering, normalization, aggregation, or reordering.

## Primary differential-accessibility analysis

Regions were retained if their count was at least 10 in at least 5 of the 10
samples. This retained 59,186 of 112,759 regions. DESeq2 1.50.2 under R 4.5.2
used the design:

```r
~ donor + cell_type
```

Mature NK was the reference level, so positive log2 fold changes represent
greater accessibility in Memory NK. The primary threshold required adjusted
p-value < 0.05 and absolute log2 fold change >= 1. Classification was performed
with ordinary DESeq2 effects and separately with apeglm-shrunken effects.

No region passed the primary threshold under either effect estimate.

VST with `blind = FALSE`, PCA, mean-SD plots, and the sample-distance heatmap
were used as exploratory diagnostics, not significance tests.

## Sensitivity and exploratory analyses

Separately labelled 10-in-3 and no-explicit-prefilter analyses retained the
same donor-aware design and contrast. All filter settings returned zero primary-
threshold DARs. A saved-table sensitivity analysis counted raw-p-value candidate
sets but did not refit DESeq2 or convert nominal candidates into discoveries.

The strict exploratory set contained 21 primary-filter regions with raw p < 0.05
and absolute apeglm LFC >= 1. GENCODE Release 50 assigned genomic categories and
nearest-gene labels. These labels are hypotheses.

g:Profiler used a custom background of 29,861 genes derived from the primary
tested peak universe. HOMER used 1,562 nominal target regions, central 200-bp
windows, and all 59,186 primary tested peaks as background. IGV reviewed six
nominal loci across all ten BAM/BAI pairs with independently autoscaled coverage
tracks. Enrichment, motif labels, and IGV views are exploratory only.
