# 10. Annotation and enrichment

[← Previous: Differential accessibility](09_differential_accessibility.md) · [Tutorial contents](README.md) · [Next: Interpretation and reporting →](11_interpretation_and_reporting.md)

## Annotation

Use an annotation release compatible with the reference assembly and chromosome
names. Promoter, intronic, exonic, and distal categories depend on declared
rules. A nearest-gene or nearest-TSS label is a proximity hypothesis, not proof
that a peak regulates that gene.

## Create coordinate-safe target and background BED files

This R block starts from the outputs of Chapters 8 and 9. It converts the
one-based inclusive SAF-style start back to a zero-based BED start. By default,
it uses only threshold-defined significant regions. If there are none, stop and
report that result unless an explicitly approved nominal exploration is being
created under a different, clearly labelled definition.

```r
results <- read.delim(
  "results/differential_accessibility/all_tested_regions.tsv",
  check.names = FALSE
)
annotation <- read.delim(
  "results/counts/consensus_peak_annotation.tsv",
  check.names = FALSE
)

colnames(annotation)[colnames(annotation) == "Geneid"] <- "peak_id"
annotation <- annotation[match(results$peak_id, annotation$peak_id), ]
if (!identical(annotation$peak_id, results$peak_id)) {
  stop("Peak annotation and result IDs do not match in order.")
}

TARGET_COLUMN <- "significant_apeglm"  # CHANGE_ME only if another rule was prespecified
if (!TARGET_COLUMN %in% colnames(results)) stop("Requested target column is absent.")
targets <- results[[TARGET_COLUMN]]
if (!any(targets)) {
  stop("No FDR-significant targets. Report the primary result before any nominal exploration.")
}

to_bed <- function(rows) {
  data.frame(
    chromosome = annotation$Chr[rows],
    start = annotation$Start[rows] - 1L,
    end = annotation$End[rows],
    peak_id = annotation$peak_id[rows]
  )
}

dir.create("results/downstream", recursive = TRUE, showWarnings = FALSE)
write.table(
  to_bed(targets),
  "results/downstream/target_regions.bed",
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
write.table(
  to_bed(rep(TRUE, nrow(results))),
  "results/downstream/tested_background_regions.bed",
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
```

## Copy-and-adapt HOMER annotation and motif commands

Replace the reference paths with resources matching the exact assembly used for
alignment. The promoter definition and motif window remain explicit choices.

```bash
set -euo pipefail

GENOME_FASTA="/path/to/reference.fa"
ANNOTATION_GTF="/path/to/build-compatible.annotation.gtf"
TARGET_BED="results/downstream/target_regions.bed"
BACKGROUND_BED="results/downstream/tested_background_regions.bed"

test -s "${GENOME_FASTA}"
test -s "${ANNOTATION_GTF}"
test -s "${TARGET_BED}"
test -s "${BACKGROUND_BED}"
mkdir -p results/downstream/annotation results/downstream/motifs logs/downstream

annotatePeaks.pl \
  "${TARGET_BED}" \
  "${GENOME_FASTA}" \
  -gtf "${ANNOTATION_GTF}" \
  > results/downstream/annotation/target_regions.homer_annotation.tsv \
  2> logs/downstream/homer_annotation.log

findMotifsGenome.pl \
  "${TARGET_BED}" \
  "${GENOME_FASTA}" \
  results/downstream/motifs \
  -size 200 \
  -bg "${BACKGROUND_BED}" \
  -p 2 \
  2>&1 | tee logs/downstream/homer_motifs.log
```

The central 200-bp window is an explicit example, not the same thing as a
promoter window. Review it for the biological question and peak widths.
Consult HOMER's official
[genomic motif-analysis documentation](https://homer.ucsd.edu/homer/ngs/peakMotifs.html)
when changing the target, background, genome, or window options.

## Enrichment background

The background should represent regions or genes that could have been selected
by the analysis. Using the whole genome when only a filtered subset was tested
can bias enrichment. Record peak-to-gene mapping, duplicate-gene handling,
database version, correction method, and the exact target and background sizes.

After producing reviewed, unique gene lists from the target annotations and
the full tested background, g:Profiler can be run from R as follows:

```r
library(gprofiler2)

TARGET_GENES_FILE <- "results/downstream/target_genes.txt"
BACKGROUND_GENES_FILE <- "results/downstream/tested_background_genes.txt"
ORGANISM_CODE <- "hsapiens"  # CHANGE_ME for another supported organism

target_genes <- unique(scan(TARGET_GENES_FILE, what = character(), quiet = TRUE))
background_genes <- unique(scan(BACKGROUND_GENES_FILE, what = character(), quiet = TRUE))

if (!length(target_genes) || !length(background_genes)) {
  stop("Target or background gene list is empty.")
}

enrichment <- gost(
  query = target_genes,
  organism = ORGANISM_CODE,
  custom_bg = background_genes,
  correction_method = "fdr",
  evcodes = TRUE
)

if (is.null(enrichment$result)) {
  message("No enrichment term passed the selected correction.")
} else {
  write.table(
    enrichment$result,
    "results/downstream/gprofiler_results.tsv",
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}
```

The g:Profiler call requires internet access. Record the query date, organism
code, source databases, correction method, target size, and background size
because the service's underlying database builds are updated over time.

## Motif analysis

Use an explicit target set and matched accessible-region background. Motif
enrichment indicates sequence-pattern enrichment; it does not establish that a
transcription factor is expressed, bound, active, or causal.

## When there are no FDR-significant regions

Lead with the zero-FDR primary result. Optional raw-p-value or effect-size
candidate sets must be labelled nominal and exploratory everywhere—in filenames,
figures, captions, tables, and prose. They must not replace the primary result.

## Downstream gate

- target definition is frozen and reported;
- background derives from the tested universe;
- genome and annotation releases match;
- multiple testing is handled within enrichment tools;
- proximity, enrichment, and motif findings use non-causal language.

---

[← Previous: Differential accessibility](09_differential_accessibility.md) · [Tutorial contents](README.md) · [Next: Interpretation and reporting →](11_interpretation_and_reporting.md)
