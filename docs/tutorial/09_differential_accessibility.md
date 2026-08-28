# 09. Differential accessibility

## Validate inputs before fitting

Confirm that matrix row identifiers are unique, matrix columns exactly equal
metadata sample IDs in the same order, counts are nonnegative integers, factor
levels are intentional, and the proposed model matrix is full rank.

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
metadata$condition <- relevel(factor(metadata$condition), ref = "condition_A")
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~ batch + condition
)
dds <- DESeq(dds, parallel = FALSE)
res <- results(
  dds,
  contrast = c("condition", "condition_B", "condition_A"),
  alpha = 0.05
)
shrunk <- lfcShrink(dds, coef = "condition_condition_B_vs_condition_A", type = "apeglm")
```

The formula, coefficient name, contrast, and reference level must be derived
from the actual metadata. Verify that the named coefficient and explicit
contrast point in the same direction.

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
