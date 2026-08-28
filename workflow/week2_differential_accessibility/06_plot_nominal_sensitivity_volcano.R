#!/usr/bin/env Rscript

options(warn = 1, scipen = 999, stringsAsFactors = FALSE, Ncpus = 1L)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLIS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

stopf <- function(...) stop(sprintf(...), call. = FALSE)
write_tsv <- function(data, path) {
  write.table(data, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_argument) != 1L) stopf("Could not resolve the script path.")
script_path <- normalizePath(sub("^--file=", "", script_argument), winslash = "/", mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)

addon_folder <- Sys.getenv(
  "VOLCANO_ADDON_FOLDER",
  unset = "volcano_multipanel_addon_attempt2"
)
addon_root <- file.path(
  project_root,
  "results", "week2_differential_accessibility",
  "nominal_sensitivity", "attempt1", "figures",
  addon_folder
)
provenance_dir <- file.path(addon_root, "provenance")
if (dir.exists(addon_root) || file.exists(addon_root)) {
  stopf("Refusing to reuse addon path: %s", addon_root)
}
dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)

write_status <- function(status, message) {
  write_tsv(
    data.frame(
      status = status,
      message = message,
      completed_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      stringsAsFactors = FALSE
    ),
    file.path(provenance_dir, "analysis_status.tsv")
  )
}

personal_library <- Sys.getenv("R_LIBS_USER")
if (!nzchar(personal_library) || !dir.exists(personal_library)) {
  write_status("FAIL", "The approved personal R library is unavailable.")
  stopf("The approved personal R library is unavailable.")
}
.libPaths(unique(c(personal_library, .libPaths())))
required_packages <- c("ggplot2", "openssl", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  message <- sprintf("Required packages are unavailable: %s", paste(missing_packages, collapse = ", "))
  write_status("FAIL", message)
  stopf("%s", message)
}

sha256_file <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  tolower(paste(format(openssl::sha256(connection)), collapse = ""))
}
relative_to_project <- function(path) {
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  prefix <- paste0(project_root, "/")
  if (startsWith(normalized, prefix)) substring(normalized, nchar(prefix) + 1L) else normalized
}
format_count <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)

scenarios <- data.frame(
  scenario_id = c("minimum_samples_5", "minimum_samples_3", "minimum_samples_0"),
  scenario_order = 1:3,
  minimum_samples = c(5L, 3L, 0L),
  scenario_plot_label = c("Minimum samples = 5\n(primary)", "Minimum samples = 3", "Minimum samples = 0\n(no explicit prefilter)"),
  expected_tested_peaks = c(59186L, 73503L, 112759L),
  expected_raw_p = c(3426L, 4079L, 4366L),
  expected_ordinary = c(1562L, 1972L, 2289L),
  expected_apeglm = c(21L, 21L, 26L),
  relative_path = c(
    "results/week2_differential_accessibility/attempt3/tables/deseq2_results_apeglm_shrunken.tsv",
    "results/week2_differential_accessibility/filter_sensitivity/10_in_3/attempt1/tables/deseq2_results_apeglm_shrunken.tsv",
    "results/week2_differential_accessibility/filter_sensitivity/no_prefilter/attempt1/tables/deseq2_results_apeglm_shrunken.tsv"
  ),
  expected_sha256 = c(
    "998363e653e10294c7c4fd2e4ca937d2a2a93b89cfad0e8f03eb654a267335d1",
    "2e83e8de28035d2a673e5b22ec9da0a634bcb30b2c7122991e2b605624c2c80f",
    "99370f331ca2d37308a12b2b6bcc922cf5ae3c5effa70118ff54b724b761acd9"
  ),
  stringsAsFactors = FALSE
)

effect_levels <- c("Ordinary LFC", "apeglm-shrunken LFC")
status_levels <- c(
  "Raw p >= 0.05",
  "Raw p < 0.05; |LFC| < 1",
  "Raw p < 0.05; |LFC| >= 1"
)
status_colors <- c(
  "Raw p >= 0.05" = "#B8C0CC",
  "Raw p < 0.05; |LFC| < 1" = "#D97706",
  "Raw p < 0.05; |LFC| >= 1" = "#2563EB"
)

plot_tables <- list()
count_tables <- list()
source_hashes <- list()
source_summaries <- list()
list_index <- 0L

tryCatch({
  for (i in seq_len(nrow(scenarios))) {
    scenario <- scenarios[i, , drop = FALSE]
    source_path <- file.path(project_root, scenario$relative_path)
    if (!file.exists(source_path)) stopf("Missing source table: %s", source_path)
    observed_hash <- sha256_file(source_path)
    if (!identical(observed_hash, scenario$expected_sha256)) {
      stopf("Source hash mismatch for %s", scenario$scenario_id)
    }

    source <- read.delim(source_path, check.names = FALSE)
    required_columns <- c(
      "peak_id", "pvalue", "padj",
      "log2FoldChange_unshrunken", "log2FoldChange_apeglm"
    )
    if (!all(required_columns %in% names(source))) {
      stopf("Required columns are missing from %s", scenario$scenario_id)
    }
    if (nrow(source) != scenario$expected_tested_peaks) {
      stopf("Unexpected row count for %s", scenario$scenario_id)
    }
    if (anyDuplicated(source$peak_id)) stopf("Duplicated peak IDs in %s", scenario$scenario_id)
    numeric_columns <- required_columns[-1]
    if (any(!vapply(source[numeric_columns], is.numeric, logical(1L)))) {
      stopf("Non-numeric result fields in %s", scenario$scenario_id)
    }
    if (any(!is.finite(as.matrix(source[numeric_columns])))) {
      stopf("Non-finite result values in %s", scenario$scenario_id)
    }
    if (any(source$pvalue <= 0 | source$pvalue > 1)) stopf("Invalid p-values in %s", scenario$scenario_id)
    if (any(source$padj < 0 | source$padj > 1)) stopf("Invalid adjusted p-values in %s", scenario$scenario_id)

    raw_flag <- source$pvalue < 0.05
    ordinary_flag <- raw_flag & abs(source$log2FoldChange_unshrunken) >= 1
    apeglm_flag <- raw_flag & abs(source$log2FoldChange_apeglm) >= 1
    fdr_ordinary <- source$padj < 0.05 & abs(source$log2FoldChange_unshrunken) >= 1
    fdr_apeglm <- source$padj < 0.05 & abs(source$log2FoldChange_apeglm) >= 1
    observed_counts <- c(sum(raw_flag), sum(ordinary_flag), sum(apeglm_flag))
    expected_counts <- c(scenario$expected_raw_p, scenario$expected_ordinary, scenario$expected_apeglm)
    if (!identical(as.integer(observed_counts), as.integer(expected_counts))) {
      stopf("Candidate-count mismatch for %s", scenario$scenario_id)
    }
    if (sum(fdr_ordinary) != 0L || sum(fdr_apeglm) != 0L) {
      stopf("Unexpected FDR-significant DARs in %s", scenario$scenario_id)
    }

    source_hashes[[i]] <- data.frame(
      scenario_id = scenario$scenario_id,
      source_path = scenario$relative_path,
      expected_sha256 = scenario$expected_sha256,
      observed_sha256 = observed_hash,
      hash_match = TRUE,
      stringsAsFactors = FALSE
    )
    source_summaries[[i]] <- data.frame(
      scenario_id = scenario$scenario_id,
      tested_peaks = nrow(source),
      minimum_pvalue = min(source$pvalue),
      maximum_minus_log10_raw_p = max(-log10(source$pvalue)),
      maximum_absolute_ordinary_lfc = max(abs(source$log2FoldChange_unshrunken)),
      maximum_absolute_apeglm_lfc = max(abs(source$log2FoldChange_apeglm)),
      stringsAsFactors = FALSE
    )

    for (effect_id in seq_along(effect_levels)) {
      list_index <- list_index + 1L
      effect_label <- effect_levels[effect_id]
      lfc <- if (effect_id == 1L) source$log2FoldChange_unshrunken else source$log2FoldChange_apeglm
      qualified <- raw_flag & abs(lfc) >= 1
      status <- ifelse(
        !raw_flag,
        status_levels[1],
        ifelse(qualified, status_levels[3], status_levels[2])
      )
      plot_tables[[list_index]] <- data.frame(
        scenario_plot_label = scenario$scenario_plot_label,
        effect_label = effect_label,
        log2_fold_change = lfc,
        minus_log10_raw_p = -log10(source$pvalue),
        status = status,
        status_order = match(status, status_levels),
        stringsAsFactors = FALSE
      )
      count_tables[[list_index]] <- data.frame(
        scenario_id = scenario$scenario_id,
        scenario_order = scenario$scenario_order,
        minimum_samples = scenario$minimum_samples,
        scenario_plot_label = scenario$scenario_plot_label,
        effect_type = effect_label,
        effect_order = effect_id,
        tested_peaks = nrow(source),
        raw_p_lt_0_05 = sum(raw_flag),
        raw_p_lt_0_05_and_abs_row_lfc_ge_1 = sum(qualified),
        fdr_lt_0_05_and_abs_row_lfc_ge_1 = if (effect_id == 1L) sum(fdr_ordinary) else sum(fdr_apeglm),
        stringsAsFactors = FALSE
      )
    }
    rm(source)
    invisible(gc())
  }

  plot_data <- do.call(rbind, plot_tables)
  panel_counts <- do.call(rbind, count_tables)
  source_hashes <- do.call(rbind, source_hashes)
  source_summaries <- do.call(rbind, source_summaries)
  rownames(plot_data) <- NULL
  rownames(panel_counts) <- NULL
  plot_data <- plot_data[order(plot_data$status_order), , drop = FALSE]
  plot_data$scenario_plot_label <- factor(
    plot_data$scenario_plot_label,
    levels = scenarios$scenario_plot_label
  )
  plot_data$effect_label <- factor(plot_data$effect_label, levels = effect_levels)
  plot_data$status <- factor(plot_data$status, levels = status_levels)
  panel_counts <- panel_counts[order(panel_counts$effect_order, panel_counts$scenario_order), , drop = FALSE]
  panel_counts$scenario_plot_label <- factor(panel_counts$scenario_plot_label, levels = scenarios$scenario_plot_label)
  panel_counts$effect_type <- factor(panel_counts$effect_type, levels = effect_levels)

  x_limit <- max(abs(plot_data$log2_fold_change))
  x_limit <- max(4.25, ceiling(x_limit * 4) / 4)
  y_limit <- max(plot_data$minus_log10_raw_p)
  y_limit <- ceiling((y_limit + 0.25) * 4) / 4
  annotations <- panel_counts
  annotations$effect_label <- factor(
    as.character(annotations$effect_type),
    levels = effect_levels
  )
  annotations$x <- -x_limit + 0.14
  annotations$y <- y_limit - 0.10
  annotations$annotation <- paste0(
    "Tested: ", format_count(annotations$tested_peaks),
    "\nRaw p <0.05: ", format_count(annotations$raw_p_lt_0_05),
    "\n+ |LFC| >=1: ", format_count(annotations$raw_p_lt_0_05_and_abs_row_lfc_ge_1)
  )

  ggplot2 <- asNamespace("ggplot2")
  volcano <- ggplot2$ggplot(
    plot_data,
    ggplot2$aes(x = log2_fold_change, y = minus_log10_raw_p, color = status)
  ) +
    ggplot2$geom_hline(
      yintercept = -log10(0.05), color = "#374151", linetype = "dashed", linewidth = 0.55
    ) +
    ggplot2$geom_vline(
      xintercept = c(-1, 1), color = "#374151", linetype = "dashed", linewidth = 0.55
    ) +
    ggplot2$geom_point(size = 0.34, alpha = 0.42, stroke = 0) +
    ggplot2$geom_label(
      data = annotations,
      mapping = ggplot2$aes(x = x, y = y, label = annotation),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3.45,
      lineheight = 1.06,
      linewidth = 0.25,
      label.padding = grid::unit(0.16, "lines"),
      fill = scales::alpha("white", 0.92),
      color = "#111827"
    ) +
    ggplot2$facet_grid(effect_label ~ scenario_plot_label) +
    ggplot2$scale_color_manual(values = status_colors, name = NULL, drop = FALSE) +
    ggplot2$scale_x_continuous(
      limits = c(-x_limit, x_limit),
      breaks = seq(-floor(x_limit), floor(x_limit), by = 1),
      expand = ggplot2$expansion(mult = c(0.01, 0.01))
    ) +
    ggplot2$scale_y_continuous(
      limits = c(0, y_limit),
      breaks = seq(0, floor(y_limit), by = 1),
      expand = ggplot2$expansion(mult = c(0, 0.01))
    ) +
    ggplot2$labs(
      title = "Raw-p-value volcano comparison across filter and LFC settings",
      subtitle = paste(
        "Columns change the explicit low-count filter; rows compare ordinary with apeglm-shrunken effects.",
        "No peak passed adjusted p <0.05 and |LFC| >=1."
      ),
      x = "log2 fold change (Memory NK relative to Mature NK)",
      y = "-log10(raw p-value)",
      caption = paste(
        "Blue: raw p <0.05 and |row-specific LFC| >=1; gold: raw p <0.05 but |LFC| <1; grey: raw p >=0.05.",
        "Dashed guides mark raw p =0.05 and |LFC| =1. These are nominal candidates, not FDR-significant DARs."
      )
    ) +
    ggplot2$theme_minimal(base_size = 13.5) +
    ggplot2$theme(
      plot.background = ggplot2$element_rect(fill = "white", color = NA),
      panel.background = ggplot2$element_rect(fill = "white", color = NA),
      panel.border = ggplot2$element_rect(fill = NA, color = "#D1D5DB", linewidth = 0.45),
      panel.grid.minor = ggplot2$element_blank(),
      panel.grid.major = ggplot2$element_line(color = "#E5E7EB", linewidth = 0.3),
      axis.text = ggplot2$element_text(color = "#374151", size = 10.5),
      axis.title = ggplot2$element_text(color = "#111827", face = "bold", size = 12.5),
      strip.background = ggplot2$element_rect(fill = "#F3F4F6", color = "#D1D5DB", linewidth = 0.5),
      strip.text = ggplot2$element_text(color = "#111827", face = "bold", size = 11.5),
      plot.title = ggplot2$element_text(color = "#111827", face = "bold", size = 22),
      plot.subtitle = ggplot2$element_text(color = "#4B5563", size = 12),
      plot.caption = ggplot2$element_text(color = "#4B5563", hjust = 0, size = 9.5),
      legend.position = "top",
      legend.justification = "left",
      legend.text = ggplot2$element_text(size = 10.5),
      legend.key.width = grid::unit(1.2, "lines"),
      panel.spacing = grid::unit(0.7, "lines"),
      plot.margin = ggplot2$margin(16, 22, 15, 16)
    )

  image_path <- file.path(addon_root, "raw_pvalue_volcano_multipanel_ordinary_and_apeglm.png")
  ggplot2$ggsave(
    filename = image_path,
    plot = volcano,
    width = 16,
    height = 9,
    units = "in",
    dpi = 300,
    bg = "white"
  )

  panel_counts_output <- panel_counts
  panel_counts_output$scenario_plot_label <- gsub(
    "\n",
    " / ",
    as.character(panel_counts_output$scenario_plot_label),
    fixed = TRUE
  )
  panel_counts_output$effect_type <- as.character(panel_counts_output$effect_type)
  write_tsv(panel_counts_output, file.path(addon_root, "volcano_panel_counts.tsv"))
  write_tsv(source_hashes, file.path(provenance_dir, "source_hashes.tsv"))
  write_tsv(source_summaries, file.path(provenance_dir, "source_summaries.tsv"))
  write_tsv(
    data.frame(
      parameter = c(
        "analysis_action", "plot_layout", "x_scale_policy", "y_scale_policy",
        "raw_p_threshold", "absolute_lfc_threshold", "fdr_threshold",
        "contrast", "minimum_samples_zero_meaning", "image_dimensions_pixels"
      ),
      value = c(
        "presentation-only plot from saved validated result tables; no DESeq2 or apeglm rerun",
        "2 rows (ordinary, apeglm) by 3 columns (minimum samples 5, 3, 0)",
        sprintf("fixed across all panels: %.2f to %.2f", -x_limit, x_limit),
        sprintf("fixed across all panels: 0 to %.2f", y_limit),
        "0.05 (nominal; not multiple-testing adjusted)",
        "1.0, applied to the LFC displayed in each row",
        "0.05; all panels have zero FDR-significant DARs at adjusted p <0.05 plus |LFC| >=1",
        "Memory NK relative to Mature NK; positive LFC means greater accessibility in Memory NK",
        "no explicit low-count prefilter; DESeq2 independent filtering remained enabled in the validated source analysis",
        "4800 x 2700"
      ),
      stringsAsFactors = FALSE
    ),
    file.path(provenance_dir, "plot_parameters.tsv")
  )
  write_tsv(
    data.frame(
      chart_id = "raw_pvalue_volcano_multipanel_ordinary_and_apeglm",
      audience = "exploratory analysis",
      question = "How do volcano patterns and nominal candidate counts differ across minimum-sample filters and ordinary versus apeglm LFCs?",
      chart_type = "six-panel volcano plot with fixed scales",
      data_grain = "one point per tested peak per panel",
      encodings = "x=row-specific LFC; y=-log10(raw p); color=nominal candidate category; facets=filter and LFC type",
      comparison_logic = "common axes and thresholds across all six panels",
      non_color_redundancy = "facets, threshold guide lines, and direct count annotations",
      stringsAsFactors = FALSE
    ),
    file.path(provenance_dir, "chart_contract.tsv")
  )
  write_tsv(
    data.frame(
      package = c("R", required_packages),
      version = c(
        as.character(getRversion()),
        vapply(required_packages, function(package) as.character(packageVersion(package)), character(1L))
      ),
      stringsAsFactors = FALSE
    ),
    file.path(provenance_dir, "package_versions.tsv")
  )
  write_tsv(
    data.frame(
      figure = basename(image_path),
      qa_status = "PENDING_MANUAL_REVIEW",
      reviewed_at = "NA",
      notes = "Inspect original-resolution PNG for clipping, overlap, legibility, honest common scales, and caveat visibility.",
      stringsAsFactors = FALSE
    ),
    file.path(provenance_dir, "manual_visual_qa_registry.tsv")
  )

  readme <- c(
    "# Six-panel raw-p-value volcano comparison",
    "",
    "This is a visual summary of the exploratory raw-p-value/apeglm sensitivity analysis.",
    "It reuses three protected result tables and does not rerun DESeq2 or apeglm.",
    "",
    "The columns show minimum_samples = 5 (primary), 3, and 0 (no explicit low-count prefilter).",
    "The rows show ordinary and apeglm-shrunken log2 fold changes. All six panels use common axes.",
    "",
    "Blue points meet raw p <0.05 and absolute row-specific LFC >=1. Gold points have raw p <0.05 but absolute LFC <1.",
    "These are nominal candidates only: all three underlying analyses retained zero DARs at adjusted p <0.05 plus absolute LFC >=1.",
    "",
    "Positive log2 fold change means greater accessibility in Memory NK relative to Mature NK."
  )
  writeLines(readme, file.path(addon_root, "README.md"), useBytes = TRUE)

  output_files <- list.files(addon_root, recursive = TRUE, full.names = TRUE)
  output_files <- output_files[
    !grepl("output_manifest\\.tsv$", output_files) &
      !grepl("analysis_status\\.tsv$", output_files) &
      !grepl("manual_visual_qa_registry\\.tsv$", output_files)
  ]
  output_manifest <- data.frame(
    relative_path = vapply(output_files, relative_to_project, character(1L)),
    bytes = as.numeric(file.info(output_files)$size),
    sha256 = vapply(output_files, sha256_file, character(1L)),
    stringsAsFactors = FALSE
  )
  output_manifest <- output_manifest[order(output_manifest$relative_path), , drop = FALSE]
  write_tsv(output_manifest, file.path(provenance_dir, "output_manifest.tsv"))
  write_status("PASS", "Six-panel volcano image created from validated saved result tables; manual visual QA remains required.")
  message("VOLCANO_MULTIPANEL_BUILD_PASS")
}, error = function(e) {
  write_status("FAIL", conditionMessage(e))
  stop(e)
})
