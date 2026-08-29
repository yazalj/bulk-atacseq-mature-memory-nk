#!/usr/bin/env Rscript

# =============================================================================
# Week 2 bulk ATAC-seq differential accessibility analysis
# Mature human NK cells versus Memory human NK cells
# =============================================================================
#
# This single script reproduces the validated Week 2 analysis from the clean
# count matrix and ordered metadata.
#
# How to run it in RStudio
# ------------------------
# 1. Open RStudio.
# 2. Set the working directory to the project root:
#    P02_NK_ATACseq_internship
# 3. Open this script and click "Source".
# 4. Review the plots in the lower-right Plots pane. Use its Back and Forward
#    arrows to move through the plots after the script finishes.
#
# Important analysis choices
# --------------------------
# - Low-count filter: retain a peak when count >= 10 in at least 5 of the
#   10 samples.
# - Design: ~ donor + cell_type
# - Reference level: Mature_NK
# - Tested contrast: Memory_NK relative to Mature_NK
# - Significant DAR: adjusted p-value < 0.05 and absolute shrunken
#   log2 fold change >= 1
# - LFC shrinkage: apeglm
# - VST: blind = FALSE
#
# Positive log2 fold changes mean greater accessibility in Memory_NK.
# Negative log2 fold changes mean greater accessibility in Mature_NK.
#
# The script never overwrites the validated analysis. It writes to a separate
# public_reproduction_run directory and stops if that directory already exists.

options(
  warn = 1,
  scipen = 999,
  stringsAsFactors = FALSE,
  Ncpus = 2L,
  repos = c(CRAN = "https://cloud.r-project.org")
)
Sys.setenv(MAKEFLAGS = "-j2")

# =============================================================================
# 1. Install and load the required packages
# =============================================================================

# Packages are installed only when they are missing. On Windows, binary
# packages are requested to avoid unnecessary compilation.
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages(
    "BiocManager",
    dependencies = NA,
    Ncpus = 2L,
    type = "binary"
  )
}

cran_packages <- c(
  "ggplot2",
  "pheatmap",
  "ggrepel",
  "hexbin",
  "openssl",
  "scales"
)

missing_cran <- cran_packages[
  !vapply(
    cran_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1L)
  )
]

if (length(missing_cran) > 0L) {
  install.packages(
    missing_cran,
    dependencies = NA,
    Ncpus = 2L,
    type = "binary"
  )
}

bioconductor_packages <- c(
  "DESeq2",
  "SummarizedExperiment",
  "apeglm",
  "vsn"
)

missing_bioconductor <- bioconductor_packages[
  !vapply(
    bioconductor_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1L)
  )
]

if (length(missing_bioconductor) > 0L) {
  BiocManager::install(
    missing_bioconductor,
    ask = FALSE,
    update = FALSE,
    dependencies = NA,
    Ncpus = 2L,
    type = "binary"
  )
}

required_packages <- c(cran_packages, bioconductor_packages)
package_available <- vapply(
  required_packages,
  requireNamespace,
  quietly = TRUE,
  FUN.VALUE = logical(1L)
)

if (!all(package_available)) {
  stop(
    "Required packages are unavailable: ",
    paste(names(package_available)[!package_available], collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(DESeq2)
  library(SummarizedExperiment)
  library(apeglm)
  library(vsn)
  library(ggplot2)
  library(pheatmap)
  library(ggrepel)
  library(hexbin)
  library(openssl)
  library(scales)
})

# =============================================================================
# 2. Project paths and small helper functions
# =============================================================================

# Keep this portable across working directories: the project root is the
# RStudio working directory rather than a personal absolute Windows path.
project_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

input_files <- c(
  count_matrix = file.path(
    project_root,
    "results",
    "week2_consensus_counts",
    "attempt1",
    "counts",
    "clean",
    "mature_memory_counts.tsv"
  ),
  peak_annotation = file.path(
    project_root,
    "results",
    "week2_consensus_counts",
    "attempt1",
    "counts",
    "clean",
    "mature_memory_peak_annotation.tsv"
  ),
  ordered_metadata = file.path(
    project_root,
    "results",
    "week2_consensus_counts",
    "attempt1",
    "counts",
    "metadata",
    "ordered_samples.tsv"
  )
)

missing_inputs <- input_files[!file.exists(input_files)]
if (length(missing_inputs) > 0L) {
  stop(
    paste0(
      "Protected inputs were not found. In RStudio, set the working ",
      "directory to the P02_NK_ATACseq_internship project root.\nMissing: ",
      paste(missing_inputs, collapse = "\n")
    ),
    call. = FALSE
  )
}

output_root <- file.path(
  project_root,
  "results",
  "week2_differential_accessibility",
  "public_reproduction_run"
)

if (dir.exists(output_root)) {
  stop(
    paste0(
      "Refusing to overwrite an earlier public reproduction run:\n",
      output_root,
      "\nRename that directory before running the script again."
    ),
    call. = FALSE
  )
}

table_dir <- file.path(output_root, "tables")
figure_dir <- file.path(output_root, "figures")
provenance_dir <- file.path(output_root, "provenance")

dir.create(table_dir, recursive = TRUE)
dir.create(figure_dir)
dir.create(provenance_dir)

check_that <- function(condition, message_text) {
  if (!isTRUE(condition)) {
    stop(message_text, call. = FALSE)
  }
}

write_tsv <- function(data, path) {
  write.table(
    data,
    path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    na = "NA"
  )
}

sha256_file <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  paste(format(openssl::sha256(connection)), collapse = "")
}

show_and_save <- function(
    plot_object,
    filename_stem,
    width,
    height
) {
  # Explicit printing sends the plot to RStudio's lower-right Plots pane.
  print(plot_object)

  ggsave(
    filename = file.path(
      figure_dir,
      paste0(filename_stem, ".png")
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  ggsave(
    filename = file.path(
      figure_dir,
      paste0(filename_stem, ".pdf")
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    bg = "white"
  )
}

# =============================================================================
# 3. Validate and load the protected inputs
# =============================================================================

expected_hashes <- c(
  count_matrix =
    "963b02f560ef1df8696b39e1189c0e9a33bf8db6e4e39f2766975570a633178f",
  peak_annotation =
    "80583d7a8359f71d43e1fe2df846607f2f5fa184287ded32af9d50563b79dca3",
  ordered_metadata =
    "dcb6b6e2e71e46f9ff002dc44f277afaeee5860462d91513884a4b66529af8cd"
)

observed_hashes <- vapply(
  input_files,
  sha256_file,
  FUN.VALUE = character(1L)
)

check_that(
  identical(observed_hashes, expected_hashes),
  "A protected input checksum has changed. Analysis stopped."
)

input_hash_table <- data.frame(
  input = names(input_files),
  path = normalizePath(
    input_files,
    winslash = "/",
    mustWork = TRUE
  ),
  expected_sha256 = unname(expected_hashes),
  observed_sha256 = unname(observed_hashes),
  match = unname(expected_hashes == observed_hashes)
)

write_tsv(
  input_hash_table,
  file.path(provenance_dir, "protected_input_hashes.tsv")
)

counts_table <- read.delim(
  input_files[["count_matrix"]],
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

peak_annotation <- read.delim(
  input_files[["peak_annotation"]],
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

metadata <- read.delim(
  input_files[["ordered_metadata"]],
  sep = "\t",
  header = TRUE,
  quote = "",
  comment.char = "",
  check.names = FALSE
)

expected_samples <- c(
  "SRR7650764",
  "SRR7650766",
  "SRR7650808",
  "SRR7650846",
  "SRR7650848",
  "SRR7650883",
  "SRR7650885",
  "SRR7650911",
  "SRR7650913",
  "SRR7650922"
)

check_that(
  identical(colnames(counts_table), c("peak_id", expected_samples)),
  "The count-matrix columns or sample order have changed."
)
check_that(
  identical(as.character(metadata$sample_id), expected_samples),
  "The ordered metadata no longer matches the validated sample order."
)
check_that(
  identical(
    as.character(counts_table$peak_id),
    as.character(peak_annotation$peak_id)
  ),
  "The count matrix and peak annotation do not have identical peak IDs."
)
check_that(
  nrow(counts_table) == 112759L,
  "The input matrix does not contain the expected 112,759 peaks."
)

count_matrix <- as.matrix(
  counts_table[, expected_samples, drop = FALSE]
)

check_that(
  !anyNA(count_matrix),
  "The count matrix contains missing values."
)
check_that(
  all(is.finite(count_matrix)),
  "The count matrix contains non-finite values."
)
check_that(
  all(count_matrix >= 0),
  "The count matrix contains negative values."
)
check_that(
  all(count_matrix == floor(count_matrix)),
  "The count matrix contains non-integer values."
)

storage.mode(count_matrix) <- "integer"
rownames(count_matrix) <- as.character(counts_table$peak_id)

metadata$sample_id <- as.character(metadata$sample_id)
metadata$donor <- factor(
  as.character(metadata$donor),
  levels = c("1001", "1002", "1003", "1004", "1008", "1010")
)
metadata$cell_type <- factor(
  as.character(metadata$cell_type),
  levels = c("Mature_NK", "Memory_NK")
)
rownames(metadata) <- metadata$sample_id

check_that(
  sum(metadata$cell_type == "Mature_NK") == 4L,
  "Expected four Mature_NK samples."
)
check_that(
  sum(metadata$cell_type == "Memory_NK") == 6L,
  "Expected six Memory_NK samples."
)
check_that(
  identical(levels(metadata$cell_type)[1L], "Mature_NK"),
  "Mature_NK is not the reference level."
)

# =============================================================================
# 4. Build the SummarizedExperiment and filter low-count peaks
# =============================================================================

nk_se <- SummarizedExperiment(
  assays = list(counts = count_matrix),
  rowData = S4Vectors::DataFrame(
    peak_annotation,
    row.names = peak_annotation$peak_id
  ),
  colData = S4Vectors::DataFrame(
    metadata,
    row.names = metadata$sample_id
  )
)

minimum_count <- 10L
minimum_samples <- 5L

keep_peak <- rowSums(
  assay(nk_se, "counts") >= minimum_count
) >= minimum_samples

nk_se_filtered <- nk_se[keep_peak, ]

check_that(
  nrow(nk_se_filtered) == 59186L,
  paste0(
    "The validated filter should retain 59,186 peaks, but retained ",
    format(nrow(nk_se_filtered), big.mark = ","),
    "."
  )
)

filtering_summary <- data.frame(
  rule = "count >= 10 in at least 5 of 10 samples",
  input_peaks = nrow(nk_se),
  retained_peaks = nrow(nk_se_filtered),
  removed_peaks = sum(!keep_peak),
  retention_percent = 100 * nrow(nk_se_filtered) / nrow(nk_se)
)

write_tsv(
  filtering_summary,
  file.path(table_dir, "filtering_summary.tsv")
)

# =============================================================================
# 5. Validate the paired donor-aware model
# =============================================================================

design_matrix <- model.matrix(
  ~ donor + cell_type,
  data = as.data.frame(colData(nk_se_filtered))
)

design_rank <- qr(design_matrix)$rank

check_that(
  identical(dim(design_matrix), c(10L, 7L)),
  "The design matrix should have dimensions 10 x 7."
)
check_that(
  design_rank == ncol(design_matrix),
  "The design matrix is not full rank."
)
check_that(
  "cell_typeMemory_NK" %in% colnames(design_matrix),
  "The Memory_NK versus Mature_NK model term is missing."
)

design_matrix_table <- data.frame(
  sample_id = rownames(design_matrix),
  design_matrix,
  check.names = FALSE
)

write_tsv(
  design_matrix_table,
  file.path(table_dir, "design_matrix.tsv")
)

# =============================================================================
# 6. Run DESeq2 and extract the requested contrast
# =============================================================================

message("Constructing DESeqDataSet with design ~ donor + cell_type.")

dds <- DESeqDataSet(
  nk_se_filtered,
  design = ~ donor + cell_type
)

message("Running DESeq2 with parallel processing disabled.")

dds <- DESeq(
  dds,
  parallel = FALSE,
  quiet = FALSE
)

size_factors <- sizeFactors(dds)

check_that(
  all(is.finite(size_factors) & size_factors > 0),
  "DESeq2 did not produce finite, positive size factors."
)

coefficient_name <- "cell_type_Memory_NK_vs_Mature_NK"

check_that(
  coefficient_name %in% resultsNames(dds),
  paste0(
    "The expected coefficient is missing. Observed coefficients: ",
    paste(resultsNames(dds), collapse = ", ")
  )
)

unshrunken_result <- results(
  dds,
  contrast = c("cell_type", "Memory_NK", "Mature_NK"),
  alpha = 0.05,
  independentFiltering = TRUE
)

# Verify that the explicit contrast has the same direction as the model
# coefficient used for apeglm shrinkage.
coefficient_result <- results(
  dds,
  name = coefficient_name,
  alpha = 0.05,
  independentFiltering = TRUE
)

finite_lfc <- (
  is.finite(unshrunken_result$log2FoldChange) &
    is.finite(coefficient_result$log2FoldChange)
)

check_that(
  any(finite_lfc),
  "No finite log2 fold changes were available for contrast checking."
)

maximum_lfc_difference <- max(
  abs(
    unshrunken_result$log2FoldChange[finite_lfc] -
      coefficient_result$log2FoldChange[finite_lfc]
  )
)

check_that(
  maximum_lfc_difference <= 1e-8,
  "The explicit contrast and the named coefficient have different LFCs."
)

shrunken_result <- lfcShrink(
  dds,
  coef = coefficient_name,
  type = "apeglm",
  parallel = FALSE
)

# =============================================================================
# 7. Prepare and save the differential-accessibility result tables
# =============================================================================

unshrunken_frame <- as.data.frame(unshrunken_result)
shrunken_frame <- as.data.frame(shrunken_result)

filtered_annotation <- peak_annotation[
  match(rownames(unshrunken_frame), peak_annotation$peak_id),
  ,
  drop = FALSE
]

check_that(
  identical(
    as.character(filtered_annotation$peak_id),
    rownames(unshrunken_frame)
  ),
  "Result rows could not be matched safely to the peak annotation."
)

significant_dar <- (
  !is.na(unshrunken_frame$padj) &
    unshrunken_frame$padj < 0.05 &
    abs(shrunken_frame$log2FoldChange) >= 1
)

dar_direction <- rep(
  "Not significant",
  nrow(unshrunken_frame)
)
dar_direction[
  significant_dar & shrunken_frame$log2FoldChange > 0
] <- "More accessible in Memory_NK"
dar_direction[
  significant_dar & shrunken_frame$log2FoldChange < 0
] <- "More accessible in Mature_NK"

unshrunken_table <- data.frame(
  filtered_annotation,
  unshrunken_frame,
  check.names = FALSE
)

shrunken_table <- data.frame(
  filtered_annotation,
  baseMean = unshrunken_frame$baseMean,
  log2FoldChange_unshrunken =
    unshrunken_frame$log2FoldChange,
  log2FoldChange_apeglm =
    shrunken_frame$log2FoldChange,
  lfcSE_apeglm = shrunken_frame$lfcSE,
  stat = unshrunken_frame$stat,
  pvalue = unshrunken_frame$pvalue,
  padj = unshrunken_frame$padj,
  classification = dar_direction,
  check.names = FALSE
)

significant_dar_table <- shrunken_table[
  significant_dar,
  ,
  drop = FALSE
]

size_factor_table <- data.frame(
  sample_id = names(size_factors),
  donor = as.character(
    metadata[names(size_factors), "donor"]
  ),
  cell_type = as.character(
    metadata[names(size_factors), "cell_type"]
  ),
  size_factor = unname(size_factors)
)

write_tsv(
  unshrunken_table,
  file.path(table_dir, "deseq2_results_unshrunken.tsv")
)
write_tsv(
  shrunken_table,
  file.path(table_dir, "deseq2_results_apeglm_shrunken.tsv")
)
write_tsv(
  significant_dar_table,
  file.path(table_dir, "significant_dars.tsv")
)
write_tsv(
  size_factor_table,
  file.path(table_dir, "size_factors.tsv")
)

# =============================================================================
# 8. Variance-stabilizing transformation and exploratory QC
# =============================================================================

vst_object <- vst(
  dds,
  blind = FALSE
)
vst_matrix <- assay(vst_object)
normalized_counts <- counts(dds, normalized = TRUE)

check_that(
  identical(dim(vst_matrix), c(59186L, 10L)),
  "The VST matrix does not have the expected dimensions."
)
check_that(
  all(is.finite(vst_matrix)),
  "The VST matrix contains non-finite values."
)

pca_data <- plotPCA(
  vst_object,
  intgroup = c("cell_type", "donor"),
  returnData = TRUE
)
pca_variance <- round(
  100 * attr(pca_data, "percentVar"),
  2
)

pca_table <- data.frame(
  sample_id = rownames(pca_data),
  PC1 = pca_data$PC1,
  PC2 = pca_data$PC2,
  cell_type = as.character(pca_data$cell_type),
  donor = as.character(pca_data$donor),
  PC1_percent_variance = pca_variance[1L],
  PC2_percent_variance = pca_variance[2L]
)

sample_distance <- as.matrix(
  dist(t(vst_matrix))
)
sample_distance_table <- data.frame(
  sample_id = rownames(sample_distance),
  sample_distance,
  check.names = FALSE
)

write_tsv(
  pca_table,
  file.path(table_dir, "pca_coordinates.tsv")
)
write_tsv(
  sample_distance_table,
  file.path(table_dir, "sample_distance_matrix.tsv")
)

# =============================================================================
# 9. Create, display, and save the plots
# =============================================================================

blue <- "#2B6CB0"
light_blue <- "#BEE3F8"
dark_blue <- "#17365D"
gold <- "#C47F00"
grey <- "#6B7280"
light_grey <- "#D1D5DB"

peak_count_plot_data <- data.frame(
  stage = factor(
    c(
      "Before low-count filtering",
      "After low-count filtering"
    ),
    levels = c(
      "Before low-count filtering",
      "After low-count filtering"
    )
  ),
  peak_count = c(
    nrow(nk_se),
    nrow(nk_se_filtered)
  )
)

peak_count_plot <- ggplot(
  peak_count_plot_data,
  aes(x = stage, y = peak_count, fill = stage)
) +
  geom_col(
    width = 0.62,
    colour = dark_blue,
    linewidth = 0.5
  ) +
  geom_text(
    aes(label = scales::comma(peak_count)),
    vjust = -0.45,
    size = 4
  ) +
  scale_fill_manual(
    values = c(light_blue, blue),
    guide = "none"
  ) +
  scale_y_continuous(
    limits = c(
      0,
      max(peak_count_plot_data$peak_count) * 1.12
    ),
    labels = scales::label_comma(),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Consensus peaks before and after filtering",
    subtitle = "Count >= 10 in at least 5 of 10 samples",
    x = NULL,
    y = "Number of consensus peaks",
    caption = paste0(
      scales::comma(sum(!keep_peak)),
      " peaks removed; ",
      round(
        100 * nrow(nk_se_filtered) / nrow(nk_se),
        2
      ),
      "% retained."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

show_and_save(
  peak_count_plot,
  "peak_count_before_after_filtering",
  width = 7,
  height = 5
)

mean_sd_before <- meanSdPlot(
  normalized_counts,
  ranks = TRUE,
  plot = FALSE,
  bins = 60
)
mean_sd_after <- meanSdPlot(
  vst_matrix,
  ranks = TRUE,
  plot = FALSE,
  bins = 60
)

make_mean_sd_plot <- function(
    mean_sd_object,
    plot_title,
    plot_subtitle,
    x_axis_title
) {
  peak_data <- data.frame(
    x = mean_sd_object$px,
    standard_deviation = mean_sd_object$py
  )
  median_data <- data.frame(
    x = mean_sd_object$rank,
    standard_deviation = mean_sd_object$sd
  )

  ggplot(
    peak_data,
    aes(x = x, y = standard_deviation)
  ) +
    geom_hex(bins = 60) +
    geom_line(
      data = median_data,
      aes(x = x, y = standard_deviation),
      inherit.aes = FALSE,
      colour = dark_blue,
      linewidth = 0.8
    ) +
    scale_fill_gradient(
      low = "#EBF8FF",
      high = blue,
      name = "Peak density"
    ) +
    scale_x_continuous(labels = scales::label_comma()) +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(
      title = plot_title,
      subtitle = plot_subtitle,
      x = x_axis_title,
      y = "Across-sample standard deviation",
      caption = paste0(
        "The dark-blue line is the running median. ",
        "All 59,186 filtered peaks are shown."
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.caption = element_text(
        colour = "#4B5563",
        margin = margin(t = 8)
      ),
      plot.margin = margin(12, 16, 22, 12)
    )
}

mean_sd_before_plot <- make_mean_sd_plot(
  mean_sd_before,
  "Mean-SD diagnostic before VST",
  paste0(
    "DESeq2 size-factor-normalized counts; ",
    "59,186 filtered peaks"
  ),
  "Rank of mean normalized count"
)

mean_sd_after_plot <- make_mean_sd_plot(
  mean_sd_after,
  "Mean-SD diagnostic after VST",
  paste0(
    "VST-transformed values with blind = FALSE; ",
    "59,186 filtered peaks"
  ),
  "Rank of mean VST value"
)

show_and_save(
  mean_sd_before_plot,
  "mean_sd_before_vst",
  width = 8.5,
  height = 6.3
)
show_and_save(
  mean_sd_after_plot,
  "mean_sd_after_vst",
  width = 8.5,
  height = 6.3
)

cell_type_colors <- c(
  Mature_NK = blue,
  Memory_NK = gold
)
donor_shapes <- c(
  `1001` = 16,
  `1002` = 17,
  `1003` = 15,
  `1004` = 18,
  `1008` = 3,
  `1010` = 4
)

pca_plot <- ggplot(
  pca_table,
  aes(
    x = PC1,
    y = PC2,
    colour = cell_type,
    shape = donor
  )
) +
  geom_hline(
    yintercept = 0,
    colour = light_grey,
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = 0,
    colour = light_grey,
    linewidth = 0.4
  ) +
  geom_point(size = 3.4, stroke = 1) +
  geom_text_repel(
    aes(label = sample_id),
    size = 3,
    colour = "#374151",
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.25,
    min.segment.length = 0,
    segment.colour = light_grey,
    seed = 20260730
  ) +
  scale_colour_manual(values = cell_type_colors) +
  scale_shape_manual(values = donor_shapes) +
  labs(
    title = "PCA of VST-transformed ATAC-seq accessibility",
    subtitle = "All 10 samples; donor-aware design; blind = FALSE",
    x = paste0(
      "PC1: ",
      format(pca_variance[1L], nsmall = 2),
      "% variance"
    ),
    y = paste0(
      "PC2: ",
      format(pca_variance[2L], nsmall = 2),
      "% variance"
    ),
    colour = "Cell type",
    shape = "Donor"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold")
  )

show_and_save(
  pca_plot,
  "pca_vst_cell_type_donor",
  width = 8.5,
  height = 6.2
)

heatmap_annotation <- data.frame(
  Cell_type = as.character(
    metadata[colnames(vst_matrix), "cell_type"]
  ),
  Donor = as.character(
    metadata[colnames(vst_matrix), "donor"]
  ),
  row.names = colnames(vst_matrix)
)

donor_colors <- c(
  `1001` = "#DEEBF7",
  `1002` = "#C6DBEF",
  `1003` = "#9ECAE1",
  `1004` = "#6BAED6",
  `1008` = "#3182BD",
  `1010` = "#08519C"
)

heatmap_object <- pheatmap(
  sample_distance,
  color = colorRampPalette(
    c("#F7FAFC", light_blue, blue, dark_blue)
  )(100),
  clustering_distance_rows = as.dist(sample_distance),
  clustering_distance_cols = as.dist(sample_distance),
  clustering_method = "complete",
  annotation_row = heatmap_annotation,
  annotation_col = heatmap_annotation,
  annotation_colors = list(
    Cell_type = cell_type_colors,
    Donor = donor_colors
  ),
  border_color = "white",
  fontsize = 9,
  fontsize_row = 8,
  fontsize_col = 8,
  angle_col = 45,
  main = paste0(
    "VST sample-to-sample Euclidean distances\n",
    "59,186 filtered peaks; blind = FALSE"
  ),
  silent = TRUE
)

# Draw the heatmap in the RStudio Plots pane.
grid::grid.newpage()
grid::grid.draw(heatmap_object$gtable)

ggsave(
  file.path(
    figure_dir,
    "sample_distance_heatmap_vst.png"
  ),
  plot = heatmap_object$gtable,
  width = 9,
  height = 8,
  units = "in",
  dpi = 300,
  bg = "white"
)
ggsave(
  file.path(
    figure_dir,
    "sample_distance_heatmap_vst.pdf"
  ),
  plot = heatmap_object$gtable,
  width = 9,
  height = 8,
  units = "in",
  bg = "white"
)

positive_padj <- unshrunken_frame$padj[
  !is.na(unshrunken_frame$padj) &
    unshrunken_frame$padj > 0
]
zero_padj_replacement <- if (length(positive_padj) > 0L) {
  max(min(positive_padj) / 10, .Machine$double.xmin)
} else {
  .Machine$double.xmin
}

plot_padj <- unshrunken_frame$padj
plot_padj[
  !is.na(plot_padj) & plot_padj == 0
] <- zero_padj_replacement

volcano_data <- data.frame(
  peak_id = filtered_annotation$peak_id,
  log2_fold_change = shrunken_frame$log2FoldChange,
  padj = unshrunken_frame$padj,
  minus_log10_padj = -log10(plot_padj),
  classification = factor(
    dar_direction,
    levels = c(
      "More accessible in Mature_NK",
      "Not significant",
      "More accessible in Memory_NK"
    )
  )
)

volcano_data <- volcano_data[
  is.finite(volcano_data$log2_fold_change) &
    is.finite(volcano_data$minus_log10_padj),
  ,
  drop = FALSE
]

volcano_labels <- volcano_data[
  volcano_data$classification != "Not significant",
  ,
  drop = FALSE
]
volcano_labels <- volcano_labels[
  order(
    volcano_labels$padj,
    -abs(volcano_labels$log2_fold_change)
  ),
  ,
  drop = FALSE
]
volcano_labels <- head(volcano_labels, 10L)

volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = log2_fold_change,
    y = minus_log10_padj,
    colour = classification
  )
) +
  geom_point(size = 1.05, alpha = 0.48) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    linewidth = 0.45
  ) +
  geom_text_repel(
    data = volcano_labels,
    aes(label = peak_id),
    size = 2.8,
    max.overlaps = Inf,
    show.legend = FALSE,
    seed = 20260730
  ) +
  scale_colour_manual(
    values = c(
      "More accessible in Mature_NK" = blue,
      "Not significant" = grey,
      "More accessible in Memory_NK" = gold
    ),
    drop = FALSE
  ) +
  labs(
    title = "Differential accessibility: Memory_NK vs Mature_NK",
    subtitle = paste0(
      "FDR < 0.05 and |apeglm-shrunken log2 fold change| >= 1"
    ),
    x = paste0(
      "apeglm-shrunken log2 fold change ",
      "(Memory_NK / Mature_NK)"
    ),
    y = "-log10(adjusted p-value)",
    colour = "Peak classification",
    caption = paste0(
      "Positive values indicate greater accessibility in Memory_NK; ",
      "negative values indicate greater accessibility in Mature_NK."
    )
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

show_and_save(
  volcano_plot,
  "volcano_memory_vs_mature_apeglm",
  width = 9.3,
  height = 7
)

# =============================================================================
# 10. Save a concise reproducibility record
# =============================================================================

analysis_summary <- data.frame(
  field = c(
    "analysis_status",
    "input_peaks",
    "retained_peaks",
    "removed_peaks",
    "filter_rule",
    "design",
    "reference_level",
    "contrast",
    "design_matrix_rank",
    "tested_peaks",
    "non_missing_adjusted_pvalues",
    "significant_dars",
    "minimum_raw_pvalue",
    "minimum_adjusted_pvalue",
    "pc1_percent_variance",
    "pc2_percent_variance"
  ),
  value = c(
    "PASS",
    nrow(nk_se),
    nrow(nk_se_filtered),
    sum(!keep_peak),
    "count >= 10 in at least 5 of 10 samples",
    "~ donor + cell_type",
    "Mature_NK",
    "Memory_NK relative to Mature_NK",
    design_rank,
    nrow(unshrunken_frame),
    sum(!is.na(unshrunken_frame$padj)),
    sum(significant_dar),
    min(unshrunken_frame$pvalue, na.rm = TRUE),
    min(unshrunken_frame$padj, na.rm = TRUE),
    pca_variance[1L],
    pca_variance[2L]
  )
)

package_versions <- data.frame(
  package = required_packages,
  version = vapply(
    required_packages,
    function(package_name) {
      as.character(packageVersion(package_name))
    },
    FUN.VALUE = character(1L)
  )
)

write_tsv(
  analysis_summary,
  file.path(provenance_dir, "analysis_summary.tsv")
)
write_tsv(
  package_versions,
  file.path(provenance_dir, "package_versions.tsv")
)
writeLines(
  capture.output(sessionInfo()),
  file.path(provenance_dir, "session_info.txt")
)

message("")
message("Analysis completed successfully.")
message("Results were written to:")
message(output_root)
message(
  "Significant DARs: ",
  nrow(significant_dar_table)
)
