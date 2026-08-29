# 09. Differential accessibility

[← Previous: Consensus peaks and counting](08_consensus_peaks_and_counting.md) · [Tutorial contents](README.md) · [Next: Annotation and enrichment →](10_annotation_and_enrichment.md)

## Validate inputs before fitting

Confirm that matrix row identifiers are unique, matrix columns exactly equal
metadata sample IDs in the same order, counts are nonnegative integers, factor
levels are intentional, and the proposed model matrix is full rank.

Start R from the project root and paste:

```r
library(DESeq2)
library(apeglm)
library(ggplot2)
library(pheatmap)

count_table <- read.delim(
  "results/counts/consensus_counts.clean.tsv",
  check.names = FALSE
)
metadata <- read.delim(
  "config/local/samples.tsv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (anyDuplicated(count_table$peak_id)) stop("Duplicate peak IDs.")
count_matrix <- as.matrix(count_table[, -1, drop = FALSE])
rownames(count_matrix) <- count_table$peak_id

if (!identical(colnames(count_matrix), metadata$sample_id)) {
  stop("Count-matrix columns and metadata sample IDs are not identical and ordered.")
}
rownames(metadata) <- metadata$sample_id
if (anyNA(count_matrix) || any(count_matrix < 0) || any(count_matrix %% 1 != 0)) {
  stop("The count matrix must contain nonnegative integers without missing values.")
}
storage.mode(count_matrix) <- "integer"
```

## Low-count filtering

Choose the rule from sample size and replicate structure before examining final
p-values. A general expression is:

```r
keep <- rowSums(count_matrix >= minimum_count) >= minimum_samples
```

The internship's `minimum_count = 10` and `minimum_samples = 5` matched its ten-
sample design. A new project must justify its own values.

## DESeq2 model pattern

```r
REFERENCE_LEVEL <- "condition_A"
NUMERATOR <- "condition_B"
DENOMINATOR <- "condition_A"
MINIMUM_COUNT <- 10L
MINIMUM_SAMPLES <- 2L  # CHANGE_ME from the prespecified replicate structure
ALPHA <- 0.05
ABSOLUTE_LFC <- 1

metadata$donor <- factor(metadata$donor)
metadata$condition <- relevel(factor(metadata$condition), ref = REFERENCE_LEVEL)

design_formula <- ~ donor + condition  # CHANGE_ME for the actual study
model_matrix <- model.matrix(design_formula, data = metadata)
if (qr(model_matrix)$rank != ncol(model_matrix)) {
  stop("Design matrix is not full rank.")
}

keep <- rowSums(count_matrix >= MINIMUM_COUNT) >= MINIMUM_SAMPLES
if (!any(keep)) stop("The low-count rule retained no regions.")
filtered_counts <- count_matrix[keep, , drop = FALSE]

dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = metadata,
  design = design_formula
)
dds <- DESeq(dds, parallel = FALSE)

res <- results(
  dds,
  contrast = c("condition", NUMERATOR, DENOMINATOR),
  alpha = ALPHA
)

print(resultsNames(dds))
COEFFICIENT_NAME <- "condition_condition_B_vs_condition_A"
if (!COEFFICIENT_NAME %in% resultsNames(dds)) {
  stop("Change COEFFICIENT_NAME to the matching value printed by resultsNames(dds).")
}

shrunk <- lfcShrink(dds, coef = COEFFICIENT_NAME, type = "apeglm")

results_table <- data.frame(
  peak_id = rownames(res),
  baseMean = res$baseMean,
  ordinary_log2FoldChange = res$log2FoldChange,
  ordinary_lfcSE = res$lfcSE,
  stat = res$stat,
  pvalue = res$pvalue,
  padj = res$padj,
  apeglm_log2FoldChange = shrunk$log2FoldChange,
  apeglm_lfcSE = shrunk$lfcSE,
  stringsAsFactors = FALSE
)

results_table$significant_ordinary <- with(
  results_table,
  !is.na(padj) & padj < ALPHA & abs(ordinary_log2FoldChange) >= ABSOLUTE_LFC
)
results_table$significant_apeglm <- with(
  results_table,
  !is.na(padj) & padj < ALPHA & abs(apeglm_log2FoldChange) >= ABSOLUTE_LFC
)

dir.create("results/differential_accessibility", recursive = TRUE, showWarnings = FALSE)
write.table(
  results_table,
  "results/differential_accessibility/all_tested_regions.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Consensus regions before filtering:", nrow(count_matrix), "\n")
cat("Regions tested:", nrow(filtered_counts), "\n")
cat("Nonmissing adjusted p-values:", sum(!is.na(res$padj)), "\n")
cat("Ordinary-LFC threshold-defined DARs:", sum(results_table$significant_ordinary), "\n")
cat("apeglm-LFC threshold-defined DARs:", sum(results_table$significant_apeglm), "\n")
```

The formula, coefficient name, contrast, and reference level must be derived
from the actual metadata. Verify that the named coefficient and explicit
contrast point in the same direction.

Create the principal exploratory QC figures from a VST object and a volcano
plot from the model results:

```r
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(
  vsd,
  intgroup = c("condition", "donor"),
  returnData = TRUE
)
percent_variance <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = donor, shape = condition)
) +
  geom_point(size = 3) +
  xlab(paste0("PC1: ", percent_variance[1], "% variance")) +
  ylab(paste0("PC2: ", percent_variance[2], "% variance")) +
  theme_bw()

ggsave(
  "results/differential_accessibility/pca_vst.png",
  pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)

sample_distances <- as.matrix(dist(t(assay(vsd))))
annotation <- metadata[, c("condition", "donor"), drop = FALSE]
rownames(annotation) <- metadata$sample_id
pheatmap(
  sample_distances,
  annotation_col = annotation,
  annotation_row = annotation,
  filename = "results/differential_accessibility/sample_distance_heatmap.png"
)

volcano_data <- results_table
volcano_data$minus_log10_padj <- -log10(
  pmax(volcano_data$padj, .Machine$double.xmin)
)

volcano_plot <- ggplot(
  volcano_data,
  aes(x = apeglm_log2FoldChange, y = minus_log10_padj, color = significant_apeglm)
) +
  geom_point(alpha = 0.6, size = 1) +
  geom_vline(xintercept = c(-ABSOLUTE_LFC, ABSOLUTE_LFC), linetype = 2) +
  geom_hline(yintercept = -log10(ALPHA), linetype = 2) +
  scale_color_manual(values = c(`FALSE` = "grey70", `TRUE` = "#0072B2")) +
  labs(
    x = paste(NUMERATOR, "relative to", DENOMINATOR, "apeglm log2 fold change"),
    y = "-log10(adjusted p-value)"
  ) +
  theme_bw()

ggsave(
  "results/differential_accessibility/volcano_apeglm.png",
  volcano_plot,
  width = 7,
  height = 5,
  dpi = 300
)
```

Review `plotDispEsts(dds)`, `sizeFactors(dds)`, and the mean-SD behavior as
additional diagnostics. A plotted separation or cluster is not a replacement
for the declared differential test.

## Review and classification

Use adjusted p-values for the primary genome-wide claim. If an effect-size
threshold is prespecified, apply it together with the adjusted-p-value rule.
Report ordinary and shrunken effects clearly. PCA, sample distances, VST plots,
size factors, and dispersion behavior are diagnostics, not tests that replace
the differential model.

## Differential-analysis gate

- design matrix is full rank;
- size factors are finite and positive;
- contrast direction is written in plain language;
- tested-region and nonmissing adjusted-p-value counts are reported;
- thresholds were fixed before classification;
- zero significant regions are reported as a valid result, not converted into
  discoveries by switching silently to raw p-values.

---

[← Previous: Consensus peaks and counting](08_consensus_peaks_and_counting.md) · [Tutorial contents](README.md) · [Next: Annotation and enrichment →](10_annotation_and_enrichment.md)
