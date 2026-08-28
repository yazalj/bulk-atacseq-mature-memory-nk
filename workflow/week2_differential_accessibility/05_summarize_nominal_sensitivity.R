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

script_argument <- grep(
  "^--file=",
  commandArgs(trailingOnly = FALSE),
  value = TRUE
)
if (length(script_argument) != 1L) {
  stopf("Could not resolve the analysis script path.")
}
script_path <- normalizePath(
  sub("^--file=", "", script_argument),
  winslash = "/",
  mustWork = TRUE
)
project_root <- normalizePath(
  file.path(dirname(script_path), "..", ".."),
  winslash = "/",
  mustWork = TRUE
)

attempt_root <- file.path(
  project_root,
  "results",
  "week2_differential_accessibility",
  "nominal_sensitivity",
  "attempt1"
)
table_dir <- file.path(attempt_root, "tables")
figure_dir <- file.path(attempt_root, "figures")
provenance_dir <- file.path(attempt_root, "provenance")

if (dir.exists(attempt_root) || file.exists(attempt_root)) {
  stopf("Refusing to reuse attempt path: %s", attempt_root)
}
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)

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
  message <- sprintf(
    "Required packages are unavailable: %s",
    paste(missing_packages, collapse = ", ")
  )
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
  if (startsWith(normalized, prefix)) {
    substring(normalized, nchar(prefix) + 1L)
  } else {
    normalized
  }
}

format_count <- function(x) format(x, big.mark = ",", scientific = FALSE)

scenarios <- data.frame(
  scenario_id = c("minimum_samples_5", "minimum_samples_3", "minimum_samples_0"),
  scenario_order = 1:3,
  minimum_samples = c(5L, 3L, 0L),
  scenario_label = c("5 (primary)", "3", "0 (no prefilter)"),
  scenario_plot_label = c("5\n(primary)", "3", "0\n(no prefilter)"),
  scenario_long_label = c(
    "Count >=10 in >=5 samples (primary)",
    "Count >=10 in >=3 samples",
    "No explicit low-count prefilter"
  ),
  filter_rule = c(
    "count >= 10 in at least 5 of 10 samples",
    "count >= 10 in at least 3 of 10 samples",
    "no explicit low-count prefilter; all protected consensus peaks supplied to DESeq2"
  ),
  expected_tested_peaks = c(59186L, 73503L, 112759L),
  apeglm_relative_path = c(
    "results/week2_differential_accessibility/attempt3/tables/deseq2_results_apeglm_shrunken.tsv",
    "results/week2_differential_accessibility/filter_sensitivity/10_in_3/attempt1/tables/deseq2_results_apeglm_shrunken.tsv",
    "results/week2_differential_accessibility/filter_sensitivity/no_prefilter/attempt1/tables/deseq2_results_apeglm_shrunken.tsv"
  ),
  unshrunken_relative_path = c(
    "results/week2_differential_accessibility/attempt3/tables/deseq2_results_unshrunken.tsv",
    "results/week2_differential_accessibility/filter_sensitivity/10_in_3/attempt1/tables/deseq2_results_unshrunken.tsv",
    "results/week2_differential_accessibility/filter_sensitivity/no_prefilter/attempt1/tables/deseq2_results_unshrunken.tsv"
  ),
  expected_apeglm_sha256 = c(
    "998363e653e10294c7c4fd2e4ca937d2a2a93b89cfad0e8f03eb654a267335d1",
    "2e83e8de28035d2a673e5b22ec9da0a634bcb30b2c7122991e2b605624c2c80f",
    "99370f331ca2d37308a12b2b6bcc922cf5ae3c5effa70118ff54b724b761acd9"
  ),
  expected_unshrunken_sha256 = c(
    "d02676374a1b58c1913254e8c9b3f0c1dec4abfdc2ecf31756742871e1f01893",
    "29566628a1f0fd41cf1a6d9875f8353923e645b662192fb4e76eee30b864c105",
    "53bf6c5aebeffd20e7addc18bafe8fef2ba638caf4060d85faa41c3572491624"
  ),
  stringsAsFactors = FALSE
)

criteria <- data.frame(
  criterion_id = c(
    "raw_p_only",
    "raw_p_and_ordinary_lfc",
    "raw_p_and_apeglm_lfc"
  ),
  criterion_order = 1:3,
  criterion_label = c(
    "Raw p <0.05",
    "Raw p <0.05 + |ordinary LFC| >=1",
    "Raw p <0.05 + |apeglm LFC| >=1"
  ),
  short_label = c(
    "Raw p only",
    "+ ordinary LFC",
    "+ apeglm LFC"
  ),
  lfc_source = c("apeglm_for_direction_only", "ordinary", "apeglm"),
  lfc_threshold = c(NA_real_, 1, 1),
  stringsAsFactors = FALSE
)

required_apeglm_columns <- c(
  "peak_id", "Chr", "Start", "End", "Strand", "Length", "baseMean",
  "log2FoldChange_unshrunken", "log2FoldChange_apeglm", "lfcSE_apeglm",
  "stat", "pvalue", "padj"
)
required_unshrunken_columns <- c(
  "peak_id", "Chr", "Start", "End", "Strand", "Length", "baseMean",
  "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
)

source_hashes <- data.frame(
  scenario_id = character(0L),
  source_role = character(0L),
  relative_path = character(0L),
  expected_sha256 = character(0L),
  observed_sha256 = character(0L),
  match = logical(0L),
  stringsAsFactors = FALSE
)
comparison_rows <- list()
scenario_rows <- list()
raw_candidate_tables <- list()
robust_candidate_tables <- list()
fdr_rows <- list()

main <- function() {
  for (i in seq_len(nrow(scenarios))) {
    scenario <- scenarios[i, , drop = FALSE]
    apeglm_path <- file.path(project_root, scenario$apeglm_relative_path)
    unshrunken_path <- file.path(project_root, scenario$unshrunken_relative_path)

    if (!file.exists(apeglm_path) || !file.exists(unshrunken_path)) {
      stopf("A required protected result table is missing for %s.", scenario$scenario_id)
    }

    observed_apeglm_hash <- sha256_file(apeglm_path)
    observed_unshrunken_hash <- sha256_file(unshrunken_path)
    source_hashes <<- rbind(
      source_hashes,
      data.frame(
        scenario_id = scenario$scenario_id,
        source_role = c("apeglm_shrunken_results", "unshrunken_results"),
        relative_path = c(
          scenario$apeglm_relative_path,
          scenario$unshrunken_relative_path
        ),
        expected_sha256 = c(
          scenario$expected_apeglm_sha256,
          scenario$expected_unshrunken_sha256
        ),
        observed_sha256 = c(observed_apeglm_hash, observed_unshrunken_hash),
        match = c(
          identical(observed_apeglm_hash, scenario$expected_apeglm_sha256),
          identical(observed_unshrunken_hash, scenario$expected_unshrunken_sha256)
        ),
        stringsAsFactors = FALSE
      )
    )
    if (!identical(observed_apeglm_hash, scenario$expected_apeglm_sha256) ||
        !identical(observed_unshrunken_hash, scenario$expected_unshrunken_sha256)) {
      stopf("Protected result hash mismatch for %s.", scenario$scenario_id)
    }

    ape <- read.delim(
      apeglm_path,
      sep = "\t",
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    ordinary <- read.delim(
      unshrunken_path,
      sep = "\t",
      header = TRUE,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )

    if (!all(required_apeglm_columns %in% names(ape))) {
      stopf("Required apeglm columns are missing for %s.", scenario$scenario_id)
    }
    if (!all(required_unshrunken_columns %in% names(ordinary))) {
      stopf("Required unshrunken columns are missing for %s.", scenario$scenario_id)
    }
    if (nrow(ape) != scenario$expected_tested_peaks ||
        nrow(ordinary) != scenario$expected_tested_peaks) {
      stopf("Unexpected tested-peak count for %s.", scenario$scenario_id)
    }
    if (anyDuplicated(ape$peak_id) || anyDuplicated(ordinary$peak_id)) {
      stopf("Duplicate peak IDs detected for %s.", scenario$scenario_id)
    }
    if (!identical(ape$peak_id, ordinary$peak_id)) {
      stopf("Peak order differs between result tables for %s.", scenario$scenario_id)
    }
    if (!identical(
      paste(ape$Chr, ape$Start, ape$End, sep = ":"),
      paste(ordinary$Chr, ordinary$Start, ordinary$End, sep = ":")
    )) {
      stopf("Peak coordinates differ between result tables for %s.", scenario$scenario_id)
    }

    numeric_columns <- c(
      "baseMean", "log2FoldChange_unshrunken", "log2FoldChange_apeglm",
      "lfcSE_apeglm", "stat", "pvalue", "padj"
    )
    if (any(!vapply(ape[numeric_columns], is.numeric, logical(1L)))) {
      stopf("Non-numeric analytical columns detected for %s.", scenario$scenario_id)
    }
    if (anyNA(ape[numeric_columns]) ||
        any(!is.finite(as.matrix(ape[numeric_columns])))) {
      stopf("Missing or non-finite analytical values detected for %s.", scenario$scenario_id)
    }
    if (any(ape$pvalue < 0 | ape$pvalue > 1) ||
        any(ape$padj < 0 | ape$padj > 1)) {
      stopf("P-values outside [0,1] detected for %s.", scenario$scenario_id)
    }

    max_differences <- c(
      pvalue = max(abs(ape$pvalue - ordinary$pvalue)),
      padj = max(abs(ape$padj - ordinary$padj)),
      stat = max(abs(ape$stat - ordinary$stat)),
      ordinary_lfc = max(abs(
        ape$log2FoldChange_unshrunken - ordinary$log2FoldChange
      ))
    )
    if (any(max_differences > 1e-12)) {
      stopf("Result-table reconciliation failed for %s.", scenario$scenario_id)
    }

    raw_p_only <- ape$pvalue < 0.05
    raw_p_ordinary <- raw_p_only & abs(ape$log2FoldChange_unshrunken) >= 1
    raw_p_apeglm <- raw_p_only & abs(ape$log2FoldChange_apeglm) >= 1
    fdr_ordinary <- ape$padj < 0.05 & abs(ape$log2FoldChange_unshrunken) >= 1
    fdr_apeglm <- ape$padj < 0.05 & abs(ape$log2FoldChange_apeglm) >= 1

    flags <- list(raw_p_only, raw_p_ordinary, raw_p_apeglm)
    direction_values <- list(
      ape$log2FoldChange_apeglm,
      ape$log2FoldChange_unshrunken,
      ape$log2FoldChange_apeglm
    )

    for (j in seq_len(nrow(criteria))) {
      selected <- flags[[j]]
      direction_lfc <- direction_values[[j]]
      comparison_rows[[length(comparison_rows) + 1L]] <<- data.frame(
        scenario_id = scenario$scenario_id,
        scenario_order = scenario$scenario_order,
        minimum_samples = scenario$minimum_samples,
        scenario_label = scenario$scenario_label,
        scenario_plot_label = scenario$scenario_plot_label,
        scenario_long_label = scenario$scenario_long_label,
        filter_rule = scenario$filter_rule,
        tested_peaks = nrow(ape),
        criterion_id = criteria$criterion_id[j],
        criterion_order = criteria$criterion_order[j],
        criterion_label = criteria$criterion_label[j],
        short_label = criteria$short_label[j],
        pvalue_threshold = 0.05,
        lfc_source = criteria$lfc_source[j],
        lfc_threshold = criteria$lfc_threshold[j],
        candidate_count = sum(selected),
        memory_more_accessible = sum(selected & direction_lfc > 0),
        mature_more_accessible = sum(selected & direction_lfc < 0),
        exactly_zero_lfc = sum(selected & direction_lfc == 0),
        percent_of_tested = 100 * sum(selected) / nrow(ape),
        stringsAsFactors = FALSE
      )
    }

    raw_candidates <- ape[raw_p_only, required_apeglm_columns, drop = FALSE]
    raw_candidates$scenario_id <- scenario$scenario_id
    raw_candidates$minimum_samples <- scenario$minimum_samples
    raw_candidates$rank_by_raw_pvalue <- rank(
      raw_candidates$pvalue,
      ties.method = "min"
    )
    raw_candidates$rank_by_absolute_apeglm_lfc <- rank(
      -abs(raw_candidates$log2FoldChange_apeglm),
      ties.method = "min"
    )
    raw_candidates$candidate_raw_p_only <- TRUE
    raw_candidates$candidate_raw_p_and_ordinary_lfc <-
      abs(raw_candidates$log2FoldChange_unshrunken) >= 1
    raw_candidates$candidate_raw_p_and_apeglm_lfc <-
      abs(raw_candidates$log2FoldChange_apeglm) >= 1
    raw_candidates$apeglm_direction <- ifelse(
      raw_candidates$log2FoldChange_apeglm > 0,
      "Memory_NK_more_accessible",
      ifelse(
        raw_candidates$log2FoldChange_apeglm < 0,
        "Mature_NK_more_accessible",
        "zero"
      )
    )
    raw_candidates <- raw_candidates[
      order(raw_candidates$pvalue, -abs(raw_candidates$log2FoldChange_apeglm)),
      ,
      drop = FALSE
    ]
    raw_candidate_tables[[scenario$scenario_id]] <<- raw_candidates

    robust_candidates <- raw_candidates[
      raw_candidates$candidate_raw_p_and_apeglm_lfc,
      ,
      drop = FALSE
    ]
    robust_candidates$rank_within_scenario <- seq_len(nrow(robust_candidates))
    robust_candidate_tables[[scenario$scenario_id]] <<- robust_candidates

    scenario_rows[[length(scenario_rows) + 1L]] <<- data.frame(
      scenario_id = scenario$scenario_id,
      minimum_samples = scenario$minimum_samples,
      filter_rule = scenario$filter_rule,
      tested_peaks = nrow(ape),
      minimum_raw_pvalue = min(ape$pvalue),
      minimum_adjusted_pvalue = min(ape$padj),
      raw_p_lt_0_05 = sum(raw_p_only),
      raw_p_lt_0_05_and_abs_ordinary_lfc_ge_1 = sum(raw_p_ordinary),
      raw_p_lt_0_05_and_abs_apeglm_lfc_ge_1 = sum(raw_p_apeglm),
      fdr_ordinary_dars = sum(fdr_ordinary),
      fdr_apeglm_dars = sum(fdr_apeglm),
      maximum_absolute_ordinary_lfc = max(abs(ape$log2FoldChange_unshrunken)),
      maximum_absolute_apeglm_lfc = max(abs(ape$log2FoldChange_apeglm)),
      max_pvalue_difference_between_tables = max_differences[["pvalue"]],
      max_padj_difference_between_tables = max_differences[["padj"]],
      max_ordinary_lfc_difference_between_tables = max_differences[["ordinary_lfc"]],
      stringsAsFactors = FALSE
    )

    fdr_rows[[length(fdr_rows) + 1L]] <<- data.frame(
      scenario_id = scenario$scenario_id,
      minimum_samples = scenario$minimum_samples,
      tested_peaks = nrow(ape),
      adjusted_pvalue_threshold = 0.05,
      absolute_lfc_threshold = 1,
      ordinary_lfc_dars = sum(fdr_ordinary),
      apeglm_lfc_dars = sum(fdr_apeglm),
      stringsAsFactors = FALSE
    )

    output_suffix <- if (scenario$minimum_samples == 0L) {
      "minimum_samples_0_no_prefilter"
    } else {
      sprintf("minimum_samples_%d", scenario$minimum_samples)
    }
    write_tsv(
      raw_candidates,
      file.path(table_dir, sprintf("raw_pvalue_candidates_%s.tsv", output_suffix))
    )
  }

  comparison <- do.call(rbind, comparison_rows)
  scenario_summary <- do.call(rbind, scenario_rows)
  fdr_summary <- do.call(rbind, fdr_rows)
  robust_all <- do.call(rbind, robust_candidate_tables)
  rownames(robust_all) <- NULL

  if (any(comparison$memory_more_accessible +
          comparison$mature_more_accessible +
          comparison$exactly_zero_lfc != comparison$candidate_count)) {
    stopf("Direction counts do not sum to candidate totals.")
  }
  for (scenario_id in scenarios$scenario_id) {
    subset_counts <- comparison$candidate_count[
      comparison$scenario_id == scenario_id
    ]
    if (!(subset_counts[1] >= subset_counts[2] &&
          subset_counts[2] >= subset_counts[3])) {
      stopf("Candidate definitions are not nested for %s.", scenario_id)
    }
  }
  if (any(fdr_summary$ordinary_lfc_dars != 0L) ||
      any(fdr_summary$apeglm_lfc_dars != 0L)) {
    stopf("A protected FDR result unexpectedly changed from zero.")
  }

  overlap_rows <- list()
  criterion_flags <- list(
    raw_p_only = "candidate_raw_p_only",
    raw_p_and_ordinary_lfc = "candidate_raw_p_and_ordinary_lfc",
    raw_p_and_apeglm_lfc = "candidate_raw_p_and_apeglm_lfc"
  )
  scenario_pairs <- combn(scenarios$scenario_id, 2L, simplify = FALSE)
  for (criterion_id in names(criterion_flags)) {
    flag_column <- criterion_flags[[criterion_id]]
    for (pair in scenario_pairs) {
      set_a <- raw_candidate_tables[[pair[1]]]
      set_b <- raw_candidate_tables[[pair[2]]]
      ids_a <- set_a$peak_id[set_a[[flag_column]]]
      ids_b <- set_b$peak_id[set_b[[flag_column]]]
      intersection_size <- length(intersect(ids_a, ids_b))
      union_size <- length(union(ids_a, ids_b))
      overlap_rows[[length(overlap_rows) + 1L]] <- data.frame(
        criterion_id = criterion_id,
        scenario_a = pair[1],
        scenario_b = pair[2],
        candidates_a = length(ids_a),
        candidates_b = length(ids_b),
        intersection = intersection_size,
        only_a = length(setdiff(ids_a, ids_b)),
        only_b = length(setdiff(ids_b, ids_a)),
        union = union_size,
        jaccard = if (union_size == 0L) NA_real_ else intersection_size / union_size,
        stringsAsFactors = FALSE
      )
    }
  }
  overlap <- do.call(rbind, overlap_rows)

  top_candidates <- do.call(
    rbind,
    lapply(raw_candidate_tables, function(x) head(x, 25L))
  )
  rownames(top_candidates) <- NULL

  comparison_export <- comparison[
    ,
    setdiff(names(comparison), "scenario_plot_label"),
    drop = FALSE
  ]
  write_tsv(
    comparison_export,
    file.path(table_dir, "candidate_definition_comparison.tsv")
  )
  write_tsv(scenario_summary, file.path(table_dir, "scenario_summary.tsv"))
  write_tsv(fdr_summary, file.path(table_dir, "fdr_result_preservation.tsv"))
  write_tsv(overlap, file.path(table_dir, "candidate_set_overlap.tsv"))
  write_tsv(
    robust_all,
    file.path(table_dir, "apeglm_robust_candidates_all_scenarios.tsv")
  )
  write_tsv(top_candidates, file.path(table_dir, "top_25_raw_pvalue_candidates_per_scenario.tsv"))
  write_tsv(source_hashes, file.path(provenance_dir, "source_file_hashes.tsv"))

  analysis_parameters <- data.frame(
    parameter = c(
      "source_analysis",
      "model",
      "contrast",
      "positive_lfc_meaning",
      "raw_pvalue_threshold",
      "ordinary_lfc_threshold",
      "apeglm_lfc_threshold",
      "filter_scenarios",
      "minimum_samples_zero_interpretation",
      "deseq2_rerun",
      "apeglm_rerun",
      "parallel_execution",
      "reporting_boundary"
    ),
    value = c(
      "validated saved Week 2 result tables",
      "~ donor + cell_type",
      "Memory_NK relative to Mature_NK",
      "greater accessibility in Memory_NK",
      "0.05",
      "absolute ordinary log2 fold change >= 1",
      "absolute apeglm-shrunken log2 fold change >= 1",
      "minimum_samples 5, 3, and 0/no explicit prefilter",
      "bypass explicit count>=10 prevalence filter; retain DESeq2 independent filtering from original fit",
      "NO",
      "NO; existing apeglm estimates reused",
      "FALSE; one R process with numerical thread limits set to 1",
      "raw-p-value results are nominal candidates, not FDR-significant DARs"
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(
    analysis_parameters,
    file.path(provenance_dir, "analysis_parameters.tsv")
  )

  package_versions <- data.frame(
    package = c("R", required_packages),
    version = c(
      as.character(getRversion()),
      vapply(required_packages, function(x) as.character(packageVersion(x)), character(1L))
    ),
    library = c(R.home(), vapply(required_packages, find.package, character(1L))),
    stringsAsFactors = FALSE
  )
  write_tsv(
    package_versions,
    file.path(provenance_dir, "package_versions.tsv")
  )

  chart_contracts <- data.frame(
    figure = c(
      "candidate_counts_by_filter_and_definition",
      "ordinary_vs_apeglm_effects_raw_p_candidates",
      "apeglm_candidate_direction_by_filter"
    ),
    analytical_question = c(
      "How do nominal candidate counts change across filter and effect-size definitions?",
      "How strongly does apeglm shrink ordinary effect estimates among raw-p-value candidates?",
      "How are the apeglm-robust candidates divided by accessibility direction?"
    ),
    supported_takeaway = c(
      "Counts depend strongly on both the explicit low-count filter and whether an ordinary or shrunken LFC threshold is required.",
      "apeglm reduces many ordinary effect estimates toward zero while preserving a smaller set with absolute shrunken LFC >=1.",
      "The final nominal candidate set contains both Memory-more and Mature-more accessible regions."
    ),
    chart_family = c("comparison", "relationship", "comparison"),
    concrete_variant = c(
      "ordered multi-series dot-line chart with log10 count axis",
      "faceted scatter with identity and effect-threshold references",
      "grouped bars with exact labels and zero baseline"
    ),
    data_rows = c(nrow(comparison), sum(comparison$candidate_count[comparison$criterion_id == "raw_p_only"]), 6L),
    renderer = "static ggplot2 PNG and PDF",
    palette_policy = c("two roots plus neutral", "one root plus neutral and gold references", "two roots plus neutral"),
    non_color_distinction = c("shape and direct labels", "point emphasis plus reference lines and facets", "direction labels and grouped position"),
    final_qa_surface = "original-resolution PNG plus PDF structural check",
    stringsAsFactors = FALSE
  )
  write_tsv(chart_contracts, file.path(provenance_dir, "chart_contracts.tsv"))

  ggplot2 <- asNamespace("ggplot2")
  comparison$scenario_plot_label <- factor(
    comparison$scenario_plot_label,
    levels = scenarios$scenario_plot_label
  )
  comparison$criterion_label <- factor(
    comparison$criterion_label,
    levels = criteria$criterion_label
  )

  palette <- c(
    "Raw p <0.05" = "#4B5563",
    "Raw p <0.05 + |ordinary LFC| >=1" = "#D97706",
    "Raw p <0.05 + |apeglm LFC| >=1" = "#2563EB"
  )
  shapes <- c(
    "Raw p <0.05" = 16,
    "Raw p <0.05 + |ordinary LFC| >=1" = 17,
    "Raw p <0.05 + |apeglm LFC| >=1" = 15
  )
  count_breaks <- c(10, 30, 100, 300, 1000, 3000, 10000)

  p_counts <- ggplot2$ggplot(
    comparison,
    ggplot2$aes(
      x = scenario_plot_label,
      y = candidate_count,
      color = criterion_label,
      shape = criterion_label,
      group = criterion_label
    )
  ) +
    ggplot2$geom_line(linewidth = 0.8) +
    ggplot2$geom_point(size = 3.4) +
    ggplot2$geom_text(
      ggplot2$aes(label = format_count(candidate_count)),
      vjust = -0.9,
      size = 3.8,
      show.legend = FALSE
    ) +
    ggplot2$scale_y_log10(
      breaks = count_breaks,
      labels = scales::label_comma()
    ) +
    ggplot2$scale_color_manual(values = palette, name = NULL) +
    ggplot2$scale_shape_manual(values = shapes, name = NULL) +
    ggplot2$labs(
      title = "Raw-p-value candidate counts across filter settings",
      subtitle = paste(
        "Donor-aware Memory NK vs Mature NK comparison; raw p <0.05;",
        "log10 count axis with exact labels"
      ),
      x = "Minimum samples with count >=10",
      y = "Candidate peaks (log10 scale)",
      caption = paste(
        "All three protected analyses retained 0 DARs at adjusted p <0.05 and |LFC| >=1.",
        "Minimum samples = 0 means no explicit low-count prefilter."
      )
    ) +
    ggplot2$theme_minimal(base_size = 13) +
    ggplot2$theme(
      plot.background = ggplot2$element_rect(fill = "white", color = NA),
      panel.background = ggplot2$element_rect(fill = "white", color = NA),
      panel.grid.minor = ggplot2$element_blank(),
      panel.grid.major.x = ggplot2$element_blank(),
      panel.grid.major.y = ggplot2$element_line(color = "#E5E7EB", linewidth = 0.35),
      axis.line.x = ggplot2$element_line(color = "#374151", linewidth = 0.45),
      axis.text = ggplot2$element_text(color = "#374151"),
      axis.title = ggplot2$element_text(color = "#111827"),
      plot.title = ggplot2$element_text(color = "#111827", face = "bold", size = 18),
      plot.subtitle = ggplot2$element_text(color = "#4B5563", size = 11.5),
      plot.caption = ggplot2$element_text(color = "#4B5563", hjust = 0, size = 9.5),
      legend.position = "top",
      legend.text = ggplot2$element_text(size = 10.5),
      plot.margin = ggplot2$margin(15, 28, 15, 15)
    )

  scatter_data <- do.call(rbind, raw_candidate_tables)
  rownames(scatter_data) <- NULL
  scatter_data$scenario_plot_label <- factor(
    scenarios$scenario_plot_label[match(scatter_data$scenario_id, scenarios$scenario_id)],
    levels = scenarios$scenario_plot_label
  )
  scatter_data$robust_status <- ifelse(
    scatter_data$candidate_raw_p_and_apeglm_lfc,
    "Raw p <0.05 + |apeglm LFC| >=1",
    "Other raw p <0.05 candidates"
  )
  scatter_limit <- max(
    abs(c(
      scatter_data$log2FoldChange_unshrunken,
      scatter_data$log2FoldChange_apeglm
    )),
    na.rm = TRUE
  )
  scatter_limit <- ceiling(scatter_limit * 2) / 2

  scatter_colors <- c(
    "Other raw p <0.05 candidates" = "#9CA3AF",
    "Raw p <0.05 + |apeglm LFC| >=1" = "#2563EB"
  )
  p_scatter <- ggplot2$ggplot(
    scatter_data,
    ggplot2$aes(
      x = log2FoldChange_unshrunken,
      y = log2FoldChange_apeglm,
      color = robust_status
    )
  ) +
    ggplot2$geom_hline(
      yintercept = c(-1, 1),
      color = "#D97706",
      linetype = "dashed",
      linewidth = 0.55
    ) +
    ggplot2$geom_vline(
      xintercept = c(-1, 1),
      color = "#D97706",
      linetype = "dashed",
      linewidth = 0.55
    ) +
    ggplot2$geom_abline(
      slope = 1,
      intercept = 0,
      color = "#374151",
      linetype = "dotted",
      linewidth = 0.7
    ) +
    ggplot2$geom_point(size = 1.15, alpha = 0.55) +
    ggplot2$facet_wrap(~scenario_plot_label, nrow = 1) +
    ggplot2$scale_color_manual(values = scatter_colors, name = NULL) +
    ggplot2$coord_equal(
      xlim = c(-scatter_limit, scatter_limit),
      ylim = c(-scatter_limit, scatter_limit),
      expand = FALSE
    ) +
    ggplot2$labs(
      title = "Ordinary and apeglm-shrunken effects among raw-p-value candidates",
      subtitle = paste(
        "Each point is a peak with raw p <0.05; gold guides mark |LFC| = 1;",
        "the dotted diagonal marks no shrinkage"
      ),
      x = "Ordinary log2 fold change",
      y = "apeglm-shrunken log2 fold change",
      caption = paste(
        "Positive values mean greater accessibility in Memory NK; negative values mean greater accessibility in Mature NK.",
        "Blue points are nominal candidates that retain |apeglm LFC| >=1."
      )
    ) +
    ggplot2$theme_minimal(base_size = 12.5) +
    ggplot2$theme(
      plot.background = ggplot2$element_rect(fill = "white", color = NA),
      panel.background = ggplot2$element_rect(fill = "white", color = NA),
      panel.grid.minor = ggplot2$element_blank(),
      panel.grid.major = ggplot2$element_line(color = "#E5E7EB", linewidth = 0.3),
      axis.line = ggplot2$element_line(color = "#374151", linewidth = 0.4),
      axis.text = ggplot2$element_text(color = "#374151"),
      axis.title = ggplot2$element_text(color = "#111827"),
      strip.text = ggplot2$element_text(face = "bold", color = "#111827"),
      plot.title = ggplot2$element_text(color = "#111827", face = "bold", size = 17),
      plot.subtitle = ggplot2$element_text(color = "#4B5563", size = 11),
      plot.caption = ggplot2$element_text(color = "#4B5563", hjust = 0, size = 9),
      legend.position = "top",
      legend.text = ggplot2$element_text(size = 10),
      plot.margin = ggplot2$margin(12, 18, 12, 12)
    )

  robust_comparison <- comparison[
    comparison$criterion_id == "raw_p_and_apeglm_lfc",
    ,
    drop = FALSE
  ]
  direction_data <- rbind(
    data.frame(
      scenario_plot_label = robust_comparison$scenario_plot_label,
      direction = "Memory NK more accessible",
      candidate_count = robust_comparison$memory_more_accessible,
      stringsAsFactors = FALSE
    ),
    data.frame(
      scenario_plot_label = robust_comparison$scenario_plot_label,
      direction = "Mature NK more accessible",
      candidate_count = robust_comparison$mature_more_accessible,
      stringsAsFactors = FALSE
    )
  )
  direction_data$scenario_plot_label <- factor(
    direction_data$scenario_plot_label,
    levels = scenarios$scenario_plot_label
  )
  direction_data$direction <- factor(
    direction_data$direction,
    levels = c("Memory NK more accessible", "Mature NK more accessible")
  )
  direction_colors <- c(
    "Memory NK more accessible" = "#2563EB",
    "Mature NK more accessible" = "#D97706"
  )
  p_direction <- ggplot2$ggplot(
    direction_data,
    ggplot2$aes(
      x = scenario_plot_label,
      y = candidate_count,
      fill = direction
    )
  ) +
    ggplot2$geom_col(
      position = ggplot2$position_dodge(width = 0.76),
      width = 0.68,
      color = "#374151",
      linewidth = 0.35
    ) +
    ggplot2$geom_text(
      ggplot2$aes(label = format_count(candidate_count)),
      position = ggplot2$position_dodge(width = 0.76),
      vjust = -0.45,
      size = 4.1
    ) +
    ggplot2$scale_fill_manual(values = direction_colors, name = NULL) +
    ggplot2$scale_y_continuous(
      expand = ggplot2$expansion(mult = c(0, 0.16)),
      labels = scales::label_comma()
    ) +
    ggplot2$labs(
      title = "Direction of raw-p-value candidates retaining |apeglm LFC| >=1",
      subtitle = "Nominal candidate counts by explicit low-count filter setting",
      x = "Minimum samples with count >=10",
      y = "Candidate peaks",
      caption = paste(
        "Selection rule: raw p <0.05 and absolute apeglm-shrunken log2 fold change >=1.",
        "These are not FDR-significant DARs."
      )
    ) +
    ggplot2$theme_minimal(base_size = 13) +
    ggplot2$theme(
      plot.background = ggplot2$element_rect(fill = "white", color = NA),
      panel.background = ggplot2$element_rect(fill = "white", color = NA),
      panel.grid.minor = ggplot2$element_blank(),
      panel.grid.major.x = ggplot2$element_blank(),
      panel.grid.major.y = ggplot2$element_line(color = "#E5E7EB", linewidth = 0.35),
      axis.line.x = ggplot2$element_line(color = "#374151", linewidth = 0.45),
      axis.text = ggplot2$element_text(color = "#374151"),
      axis.title = ggplot2$element_text(color = "#111827"),
      plot.title = ggplot2$element_text(color = "#111827", face = "bold", size = 17),
      plot.subtitle = ggplot2$element_text(color = "#4B5563", size = 11.5),
      plot.caption = ggplot2$element_text(color = "#4B5563", hjust = 0, size = 9.5),
      legend.position = "top",
      legend.text = ggplot2$element_text(size = 10.5),
      plot.margin = ggplot2$margin(14, 24, 14, 14)
    )

  save_plot_pair <- function(plot, stem, width = 13.333, height = 7.5) {
    ggplot2$ggsave(
      filename = file.path(figure_dir, paste0(stem, ".png")),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = 300,
      bg = "white"
    )
    ggplot2$ggsave(
      filename = file.path(figure_dir, paste0(stem, ".pdf")),
      plot = plot,
      width = width,
      height = height,
      units = "in",
      device = grDevices::cairo_pdf,
      bg = "white"
    )
  }

  save_plot_pair(p_counts, "candidate_counts_by_filter_and_definition")
  save_plot_pair(p_scatter, "ordinary_vs_apeglm_effects_raw_p_candidates")
  save_plot_pair(p_direction, "apeglm_candidate_direction_by_filter")

  visual_registry <- data.frame(
    figure_stem = c(
      "candidate_counts_by_filter_and_definition",
      "ordinary_vs_apeglm_effects_raw_p_candidates",
      "apeglm_candidate_direction_by_filter"
    ),
    png_path = file.path(
      "figures",
      paste0(c(
        "candidate_counts_by_filter_and_definition",
        "ordinary_vs_apeglm_effects_raw_p_candidates",
        "apeglm_candidate_direction_by_filter"
      ), ".png")
    ),
    pdf_path = file.path(
      "figures",
      paste0(c(
        "candidate_counts_by_filter_and_definition",
        "ordinary_vs_apeglm_effects_raw_p_candidates",
        "apeglm_candidate_direction_by_filter"
      ), ".pdf")
    ),
    qa_status = "PENDING_MANUAL_REVIEW",
    reviewed_at = "NA",
    notes = "Inspect original-resolution PNG for clipping, overlap, scale honesty, and legibility; structurally verify PDF.",
    stringsAsFactors = FALSE
  )
  write_tsv(
    visual_registry,
    file.path(provenance_dir, "manual_visual_qa_registry.tsv")
  )

  readme_lines <- c(
    "# Tutor-directed raw-p-value and apeglm sensitivity comparison",
    "",
    "## Outcome",
    "",
    paste0(
      "This non-overwriting comparison reuses the three validated Week 2 result tables. ",
      "It does not rerun DESeq2 or apeglm and does not replace the protected primary 10-in-5 analysis."
    ),
    "",
    "All three scenarios still contain zero DARs at adjusted p-value <0.05 plus absolute LFC >=1.",
    "The raw-p-value rows below are exploratory nominal candidate sets, not FDR-significant discoveries.",
    "",
    "## Candidate-count comparison",
    "",
    "| Minimum samples | Tested peaks | Raw p <0.05 | Raw p + abs(ordinary LFC) >=1 | Raw p + abs(apeglm LFC) >=1 |",
    "|---:|---:|---:|---:|---:|",
    vapply(seq_len(nrow(scenario_summary)), function(i) {
      paste0(
        "| ", scenario_summary$minimum_samples[i],
        " | ", format_count(scenario_summary$tested_peaks[i]),
        " | ", format_count(scenario_summary$raw_p_lt_0_05[i]),
        " | ", format_count(scenario_summary$raw_p_lt_0_05_and_abs_ordinary_lfc_ge_1[i]),
        " | ", format_count(scenario_summary$raw_p_lt_0_05_and_abs_apeglm_lfc_ge_1[i]),
        " |"
      )
    }, character(1L)),
    "",
    "Minimum samples = 0 means that the explicit count>=10 prevalence prefilter was bypassed; it does not disable DESeq2's independent filtering in the saved analysis.",
    "",
    "## Presentation boundary",
    "",
    "Use the phrase **nominal raw-p-value candidates with apeglm-stabilized effect estimates**.",
    "Do not call these peaks FDR-significant DARs, and keep the adjusted p-values visible in tables.",
    "Positive LFC means greater accessibility in Memory NK relative to Mature NK.",
    "",
    "## Main files",
    "",
    "- `tables/candidate_definition_comparison.tsv`: exact 3-by-3 comparison.",
    "- `tables/scenario_summary.tsv`: one-row summary per filter setting.",
    "- `tables/candidate_set_overlap.tsv`: identity overlap across filter settings.",
    "- `tables/apeglm_robust_candidates_all_scenarios.tsv`: combined nominal candidates passing raw p <0.05 and |apeglm LFC| >=1.",
    "- `figures/candidate_counts_by_filter_and_definition.*`: headline comparison.",
    "- `figures/ordinary_vs_apeglm_effects_raw_p_candidates.*`: shrinkage diagnostic.",
    "- `figures/apeglm_candidate_direction_by_filter.*`: direction counts."
  )
  writeLines(readme_lines, file.path(attempt_root, "README.md"), useBytes = TRUE)

  command_record <- c(
    paste("Rscript", shQuote(relative_to_project(script_path))),
    "Rerun policy: REFUSE if attempt root exists",
    "DESeq2 rerun: NO",
    "apeglm rerun: NO",
    "Thread policy: one R process; numerical thread environment variables set to 1"
  )
  writeLines(
    command_record,
    file.path(provenance_dir, "execution_command.txt"),
    useBytes = TRUE
  )
  capture.output(
    sessionInfo(),
    file = file.path(provenance_dir, "session_info.txt")
  )

  output_files <- list.files(
    attempt_root,
    recursive = TRUE,
    full.names = TRUE,
    all.files = FALSE
  )
  output_files <- output_files[file.info(output_files)$isdir %in% FALSE]
  output_files <- output_files[
    basename(output_files) != "output_sha256_manifest.tsv" &
      basename(output_files) != "analysis_status.tsv" &
      basename(output_files) != "manual_visual_qa_registry.tsv"
  ]
  output_manifest <- data.frame(
    relative_path = vapply(output_files, relative_to_project, character(1L)),
    bytes = file.info(output_files)$size,
    sha256 = vapply(output_files, sha256_file, character(1L)),
    stringsAsFactors = FALSE
  )
  output_manifest <- output_manifest[order(output_manifest$relative_path), , drop = FALSE]
  write_tsv(
    output_manifest,
    file.path(provenance_dir, "output_sha256_manifest.tsv")
  )

  write_status(
    "PASS_PENDING_MANUAL_VISUAL_QA",
    paste(
      "Sensitivity tables and figures created from validated saved results;",
      "DESeq2 and apeglm were not rerun."
    )
  )

  cat("TUTOR_RAW_PVALUE_APEGLM_SENSITIVITY_ANALYSIS=PASS_PENDING_MANUAL_VISUAL_QA\n")
  invisible(TRUE)
}

tryCatch(
  main(),
  error = function(e) {
    try(write_status("FAIL", conditionMessage(e)), silent = TRUE)
    stop(e)
  }
)
