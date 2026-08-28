#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: 01_annotate_strict_candidates_attempt2.R <project_root> <gencode_v50_gtf_gz>")
}

project_root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
gtf_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)

personal_library <- Sys.getenv("R_LIBS_USER", unset = "")
if (nzchar(personal_library) && dir.exists(personal_library)) {
  .libPaths(c(personal_library, .libPaths()))
}

required_packages <- c(
  "digest", "GenomicRanges", "IRanges", "GenomeInfoDb", "rtracklayer",
  "ggplot2"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
  stop("Required installed packages are unavailable: ", paste(missing_packages, collapse = ", "))
}

result_root <- file.path(project_root, "results", "week3_downstream", "annotation", "attempt1")
if (file.exists(result_root) || dir.exists(result_root)) {
  stop("Refusing to overwrite existing Week 3 attempt: ", result_root)
}

output_dirs <- file.path(
  result_root,
  c("bed", "tables", "figures", "provenance")
)
for (path in output_dirs) {
  if (!dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create output directory: ", path)
  }
}

bed_dir <- file.path(result_root, "bed")
table_dir <- file.path(result_root, "tables")
figure_dir <- file.path(result_root, "figures")
provenance_dir <- file.path(result_root, "provenance")

candidate_path <- file.path(
  project_root, "results", "week2_differential_accessibility",
  "nominal_sensitivity", "attempt1", "tables",
  "raw_pvalue_candidates_minimum_samples_5.tsv"
)
ordinary_path <- file.path(
  project_root, "results", "week2_differential_accessibility", "attempt3",
  "tables", "deseq2_results_unshrunken.tsv"
)
shrunken_path <- file.path(
  project_root, "results", "week2_differential_accessibility", "attempt3",
  "tables", "deseq2_results_apeglm_shrunken.tsv"
)
counts_path <- file.path(
  project_root, "results", "week2_consensus_counts", "attempt1", "counts",
  "clean", "mature_memory_counts.tsv"
)
metadata_path <- file.path(
  project_root, "results", "week2_consensus_counts", "attempt1", "counts",
  "metadata", "ordered_samples.tsv"
)
size_factor_path <- file.path(
  project_root, "results", "week2_differential_accessibility", "attempt3",
  "tables", "size_factors.tsv"
)

input_paths <- c(
  primary_candidate_table = candidate_path,
  primary_unshrunken_results = ordinary_path,
  primary_apeglm_results = shrunken_path,
  protected_counts = counts_path,
  protected_metadata = metadata_path,
  saved_size_factors = size_factor_path,
  gencode_v50_primary_assembly_gtf = gtf_path
)

missing_inputs <- input_paths[!file.exists(input_paths)]
if (length(missing_inputs) > 0L) {
  stop("Required input files are missing: ", paste(missing_inputs, collapse = ", "))
}

expected_sha256 <- c(
  primary_candidate_table = "9218ff76ccf340fa027faa1ddd9ec09969b0f86035be704cdb2838cb14b8092e",
  primary_unshrunken_results = "d02676374a1b58c1913254e8c9b3f0c1dec4abfdc2ecf31756742871e1f01893",
  primary_apeglm_results = "998363e653e10294c7c4fd2e4ca937d2a2a93b89cfad0e8f03eb654a267335d1",
  protected_counts = "963b02f560ef1df8696b39e1189c0e9a33bf8db6e4e39f2766975570a633178f",
  protected_metadata = "dcb6b6e2e71e46f9ff002dc44f277afaeee5860462d91513884a4b66529af8cd",
  saved_size_factors = "6a3144ab2b188aae40066c69d84436961ebe8c2c4ef82ca25afd70d8374594f7",
  gencode_v50_primary_assembly_gtf = "89bbad69a8c89fee5fadec0a2f14a098752d27a4ef7f24e4de9bac681e1b18f4"
)

sha256_file <- function(path) {
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}

actual_sha256 <- vapply(input_paths, sha256_file, FUN.VALUE = character(1))
hash_match <- actual_sha256 == expected_sha256[names(actual_sha256)]
input_hashes <- data.frame(
  input_id = names(input_paths),
  path = unname(input_paths),
  expected_sha256 = unname(expected_sha256[names(input_paths)]),
  actual_sha256 = unname(actual_sha256),
  match = unname(hash_match)
)
write.table(
  input_hashes,
  file.path(provenance_dir, "input_sha256_verification.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = ""
)
if (!all(hash_match)) {
  stop("One or more protected/reference input hashes do not match")
}

read_tsv <- function(path, check.names = FALSE) {
  read.delim(
    path, sep = "\t", header = TRUE, quote = "", comment.char = "",
    check.names = check.names, stringsAsFactors = FALSE
  )
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

write_bed <- function(x, path) {
  write.table(
    x, path, sep = "\t", quote = FALSE, row.names = FALSE,
    col.names = FALSE, na = ""
  )
}

as_flag <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  toupper(as.character(x)) == "TRUE"
}

chromosome_rank <- function(chr) {
  canonical <- c(as.character(1:22), "X", "Y", "MT", "M")
  idx <- match(chr, canonical)
  idx[is.na(idx)] <- 1000L
  idx
}

sort_intervals <- function(x) {
  x[order(chromosome_rank(x$Chr), x$Chr, x$Start, x$End, x$peak_id), , drop = FALSE]
}

to_gencode_seqname <- function(chr) {
  chr <- as.character(chr)
  out <- ifelse(chr %in% c("MT", "M"), "chrM", chr)
  ifelse(startsWith(out, "chr"), out, paste0("chr", out))
}

candidate_rows <- read_tsv(candidate_path)
ordinary_rows <- read_tsv(ordinary_path)
shrunken_rows <- read_tsv(shrunken_path)

if (nrow(candidate_rows) != 3426L) {
  stop("Expected 3,426 primary raw-p candidates; observed ", nrow(candidate_rows))
}
if (nrow(ordinary_rows) != 59186L || nrow(shrunken_rows) != 59186L) {
  stop("Protected primary result tables must each contain 59,186 rows")
}
if (!identical(ordinary_rows$peak_id, shrunken_rows$peak_id)) {
  stop("Ordinary and apeglm primary peak orders differ")
}

strict <- candidate_rows[as_flag(candidate_rows$candidate_raw_p_and_apeglm_lfc), , drop = FALSE]
broader <- candidate_rows[as_flag(candidate_rows$candidate_raw_p_and_ordinary_lfc), , drop = FALSE]
if (nrow(strict) != 21L || nrow(broader) != 1562L) {
  stop("Candidate definitions did not reproduce the approved 21 and 1,562 row sets")
}
if (sum(strict$apeglm_direction == "Memory_NK_more_accessible") != 11L ||
    sum(strict$apeglm_direction == "Mature_NK_more_accessible") != 10L) {
  stop("Strict candidate direction split did not reproduce 11 Memory-more and 10 Mature-more")
}
if (anyDuplicated(strict$peak_id) || anyDuplicated(broader$peak_id) ||
    anyDuplicated(ordinary_rows$peak_id)) {
  stop("Peak IDs must be unique within every target/background set")
}
if (any(strict$Start < 1L) || any(strict$End < strict$Start)) {
  stop("Strict one-based peak coordinates are invalid")
}

make_bed <- function(x) {
  out <- data.frame(
    chrom = x$Chr,
    chromStart = as.integer(x$Start) - 1L,
    chromEnd = as.integer(x$End),
    name = x$peak_id,
    score = 0L,
    strand = "."
  )
  if (any(out$chromStart < 0L) || any(out$chromEnd <= out$chromStart)) {
    stop("BED conversion produced invalid intervals")
  }
  out
}

strict_bed_source <- sort_intervals(strict[, c("peak_id", "Chr", "Start", "End")])
broader_bed_source <- sort_intervals(broader[, c("peak_id", "Chr", "Start", "End")])
background_bed_source <- sort_intervals(ordinary_rows[, c("peak_id", "Chr", "Start", "End")])
strict_bed <- make_bed(strict_bed_source)
broader_bed <- make_bed(broader_bed_source)
background_bed <- make_bed(background_bed_source)

write_bed(strict_bed, file.path(bed_dir, "strict_21_primary_nominal_candidates.bed"))
write_bed(broader_bed, file.path(bed_dir, "broader_1562_homer_sensitivity_candidates.bed"))
write_bed(background_bed, file.path(bed_dir, "primary_59186_tested_peaks_background.bed"))

strict_roundtrip <- data.frame(
  peak_id = strict_bed$name,
  original_start = strict_bed_source$Start,
  bed_start = strict_bed$chromStart,
  reconstructed_start = strict_bed$chromStart + 1L,
  original_end = strict_bed_source$End,
  bed_end = strict_bed$chromEnd,
  start_match = strict_bed_source$Start == strict_bed$chromStart + 1L,
  end_match = strict_bed_source$End == strict_bed$chromEnd
)
write_tsv(strict_roundtrip, file.path(table_dir, "strict_candidate_bed_coordinate_roundtrip.tsv"))
if (!all(strict_roundtrip$start_match) || !all(strict_roundtrip$end_match)) {
  stop("Strict BED coordinate round-trip validation failed")
}

message("Importing GENCODE Release 50 transcript, exon, UTR, and CDS features")
selected_feature_types <- c("transcript", "exon", "UTR", "CDS")
gtf <- rtracklayer::import(
  gtf_path,
  format = "gtf",
  feature.type = selected_feature_types
)
feature_type <- as.character(S4Vectors::mcols(gtf)$type)
if (!all(selected_feature_types %in% unique(feature_type))) {
  stop("GENCODE import did not contain every required feature type")
}
available_gencode_seqnames <- unique(as.character(GenomicRanges::seqnames(gtf)))

strict_seqnames <- unique(to_gencode_seqname(strict$Chr))
gtf <- gtf[as.character(GenomicRanges::seqnames(gtf)) %in% strict_seqnames]
feature_type <- as.character(S4Vectors::mcols(gtf)$type)
if (length(gtf) == 0L) {
  stop("No GENCODE features matched the strict candidate chromosomes")
}

required_attributes <- c("gene_id", "transcript_id", "gene_type", "gene_name", "transcript_type")
missing_attributes <- setdiff(required_attributes, colnames(S4Vectors::mcols(gtf)))
if (length(missing_attributes) > 0L) {
  stop("GENCODE import is missing required attributes: ", paste(missing_attributes, collapse = ", "))
}

tx <- gtf[feature_type == "transcript"]
exon <- gtf[feature_type == "exon"]
utr <- gtf[feature_type == "UTR"]
cds <- gtf[feature_type == "CDS"]

# GENCODE GTF uses a single UTR feature label. Derive five-prime and
# three-prime UTR segments from the transcript strand and CDS span.
cds_transcript_id <- as.character(S4Vectors::mcols(cds)$transcript_id)
utr_transcript_id <- as.character(S4Vectors::mcols(utr)$transcript_id)
cds_df <- data.frame(
  transcript_id = cds_transcript_id,
  cds_start = GenomicRanges::start(cds),
  cds_end = GenomicRanges::end(cds)
)
cds_min <- aggregate(cds_start ~ transcript_id, data = cds_df, min)
cds_max <- aggregate(cds_end ~ transcript_id, data = cds_df, max)
cds_span <- merge(cds_min, cds_max, by = "transcript_id", sort = FALSE)
utr_cds_idx <- match(utr_transcript_id, cds_span$transcript_id)
utr_has_cds <- !is.na(utr_cds_idx)
utr_strand <- as.character(GenomicRanges::strand(utr))
utr_start <- GenomicRanges::start(utr)
utr_end <- GenomicRanges::end(utr)
matched_cds_start <- rep(NA_integer_, length(utr))
matched_cds_end <- rep(NA_integer_, length(utr))
matched_cds_start[utr_has_cds] <- cds_span$cds_start[utr_cds_idx[utr_has_cds]]
matched_cds_end[utr_has_cds] <- cds_span$cds_end[utr_cds_idx[utr_has_cds]]
is_5utr <- utr_has_cds & (
  (utr_strand == "+" & utr_end < matched_cds_start) |
    (utr_strand == "-" & utr_start > matched_cds_end)
)
is_3utr <- utr_has_cds & (
  (utr_strand == "+" & utr_start > matched_cds_end) |
    (utr_strand == "-" & utr_end < matched_cds_start)
)
utr5 <- utr[is_5utr]
utr3 <- utr[is_3utr]
ambiguous_utr_segments <- sum(!(is_5utr | is_3utr))
if (length(utr5) == 0L || length(utr3) == 0L) {
  stop("Unable to derive both five-prime and three-prime UTR features from GENCODE GTF")
}
rm(gtf)
invisible(gc())

tx_strand <- as.character(GenomicRanges::strand(tx))
tx_tss <- ifelse(tx_strand == "+", GenomicRanges::start(tx), GenomicRanges::end(tx))
tx_df <- data.frame(
  seqname = as.character(GenomicRanges::seqnames(tx)),
  start = GenomicRanges::start(tx),
  end = GenomicRanges::end(tx),
  strand = tx_strand,
  tss = as.integer(tx_tss),
  gene_id = as.character(S4Vectors::mcols(tx)$gene_id),
  gene_id_unversioned = sub("\\..*$", "", as.character(S4Vectors::mcols(tx)$gene_id)),
  gene_name = as.character(S4Vectors::mcols(tx)$gene_name),
  gene_type = as.character(S4Vectors::mcols(tx)$gene_type),
  transcript_id = as.character(S4Vectors::mcols(tx)$transcript_id),
  transcript_type = as.character(S4Vectors::mcols(tx)$transcript_type)
)
if (nrow(tx_df) == 0L || anyDuplicated(tx_df$transcript_id)) {
  stop("Transcript feature table is empty or transcript IDs are duplicated")
}

strict_annot_order <- strict[order(strict$rank_by_raw_pvalue), , drop = FALSE]
peak_gr <- GenomicRanges::GRanges(
  seqnames = to_gencode_seqname(strict_annot_order$Chr),
  ranges = IRanges::IRanges(start = strict_annot_order$Start, end = strict_annot_order$End),
  strand = "*"
)

tss_gr <- GenomicRanges::GRanges(
  seqnames = tx_df$seqname,
  ranges = IRanges::IRanges(start = tx_df$tss, end = tx_df$tss),
  strand = tx_df$strand,
  tx_index = seq_len(nrow(tx_df))
)

promoter_start <- pmax(1L, tx_df$tss - 3000L)
promoter_end <- tx_df$tss + 3000L
promoter_gr <- GenomicRanges::GRanges(
  seqnames = tx_df$seqname,
  ranges = IRanges::IRanges(start = promoter_start, end = promoter_end),
  strand = tx_df$strand
)

downstream_start <- ifelse(tx_df$strand == "+", tx_df$end + 1L, pmax(1L, tx_df$start - 3000L))
downstream_end <- ifelse(tx_df$strand == "+", tx_df$end + 3000L, tx_df$start - 1L)
valid_downstream <- downstream_end >= downstream_start
downstream_gr <- GenomicRanges::GRanges(
  seqnames = tx_df$seqname[valid_downstream],
  ranges = IRanges::IRanges(
    start = downstream_start[valid_downstream],
    end = downstream_end[valid_downstream]
  ),
  strand = tx_df$strand[valid_downstream]
)

overlaps_any <- function(query, subject) {
  if (length(subject) == 0L) {
    return(rep(FALSE, length(query)))
  }
  GenomicRanges::countOverlaps(query, subject, ignore.strand = TRUE) > 0L
}

overlap_promoter <- overlaps_any(peak_gr, promoter_gr)
overlap_utr5 <- overlaps_any(peak_gr, utr5)
overlap_utr3 <- overlaps_any(peak_gr, utr3)
overlap_exon <- overlaps_any(peak_gr, exon)
overlap_transcript <- overlaps_any(peak_gr, tx)
overlap_downstream <- overlaps_any(peak_gr, downstream_gr)

annotation_category <- ifelse(
  overlap_promoter, "Promoter",
  ifelse(
    overlap_utr5, "5UTR",
    ifelse(
      overlap_utr3, "3UTR",
      ifelse(
        overlap_exon, "Exon",
        ifelse(
          overlap_transcript, "Intron",
          ifelse(overlap_downstream, "Downstream", "Distal/intergenic")
        )
      )
    )
  )
)

nearest_hits <- GenomicRanges::distanceToNearest(
  peak_gr, tss_gr, ignore.strand = TRUE, select = "all"
)
if (length(unique(S4Vectors::queryHits(nearest_hits))) != length(peak_gr)) {
  stop("At least one strict candidate lacks a nearest transcript TSS")
}

signed_tss_distance <- function(peak_start, peak_end, tss, strand) {
  if (tss >= peak_start && tss <= peak_end) {
    return(0L)
  }
  if (strand == "+") {
    if (peak_end < tss) peak_end - tss else peak_start - tss
  } else {
    if (peak_end < tss) tss - peak_end else tss - peak_start
  }
}

nearest_rows <- lapply(seq_len(length(peak_gr)), function(i) {
  hit_idx <- which(S4Vectors::queryHits(nearest_hits) == i)
  tx_idx <- S4Vectors::subjectHits(nearest_hits)[hit_idx]
  hit_tx <- tx_df[tx_idx, , drop = FALSE]
  signed_distance <- vapply(
    seq_len(nrow(hit_tx)),
    function(j) signed_tss_distance(
      GenomicRanges::start(peak_gr)[i], GenomicRanges::end(peak_gr)[i],
      hit_tx$tss[j], hit_tx$strand[j]
    ),
    FUN.VALUE = integer(1)
  )
  protein_priority <- hit_tx$transcript_type != "protein_coding"
  missing_symbol <- is.na(hit_tx$gene_name) | hit_tx$gene_name == ""
  ord <- order(abs(signed_distance), protein_priority, missing_symbol, hit_tx$gene_name, hit_tx$transcript_id)
  hit_tx <- hit_tx[ord, , drop = FALSE]
  signed_distance <- signed_distance[ord]
  primary <- hit_tx[1L, , drop = FALSE]
  unique_gene_ids <- unique(hit_tx$gene_id)
  unique_gene_symbols <- unique(hit_tx$gene_name[!is.na(hit_tx$gene_name) & hit_tx$gene_name != ""])
  data.frame(
    nearest_gene_id = primary$gene_id,
    nearest_gene_id_unversioned = primary$gene_id_unversioned,
    nearest_gene_symbol = primary$gene_name,
    nearest_gene_type = primary$gene_type,
    nearest_transcript_id = primary$transcript_id,
    nearest_transcript_type = primary$transcript_type,
    nearest_transcript_strand = primary$strand,
    nearest_tss = primary$tss,
    signed_distance_to_tss_bp = signed_distance[1L],
    nearest_tss_tie_transcripts = nrow(hit_tx),
    nearest_tss_tie_genes = length(unique_gene_ids),
    nearest_gene_ambiguous = length(unique_gene_ids) > 1L,
    tied_gene_ids = paste(unique_gene_ids, collapse = ";"),
    tied_gene_symbols = paste(unique_gene_symbols, collapse = ";")
  )
})
nearest_df <- do.call(rbind, nearest_rows)

annotation <- cbind(
  strict_annot_order,
  data.frame(
    annotation_category = annotation_category,
    overlaps_promoter_pm3kb = overlap_promoter,
    overlaps_5utr = overlap_utr5,
    overlaps_3utr = overlap_utr3,
    overlaps_exon = overlap_exon,
    overlaps_transcript = overlap_transcript,
    overlaps_downstream_3kb = overlap_downstream
  ),
  nearest_df
)

if (nrow(annotation) != 21L || anyDuplicated(annotation$peak_id) ||
    any(is.na(annotation$nearest_gene_id)) || any(annotation$nearest_gene_id == "")) {
  stop("Final strict annotation table failed row, uniqueness, or nearest-gene checks")
}
write_tsv(annotation, file.path(table_dir, "strict_21_candidate_annotation.tsv"))

counts <- read_tsv(counts_path)
metadata <- read_tsv(metadata_path)
size_factors <- read_tsv(size_factor_path)
if (!identical(colnames(counts)[-1L], metadata$sample_id) ||
    !identical(metadata$sample_id, size_factors$sample_id)) {
  stop("Count, metadata, and saved-size-factor sample orders do not match")
}
if (any(!is.finite(size_factors$size_factor)) || any(size_factors$size_factor <= 0)) {
  stop("Saved DESeq2 size factors must be positive and finite")
}

strict_count_idx <- match(annotation$peak_id, counts$peak_id)
if (anyNA(strict_count_idx)) {
  stop("At least one strict candidate is absent from the protected count matrix")
}
count_matrix <- as.matrix(counts[strict_count_idx, -1L, drop = FALSE])
storage.mode(count_matrix) <- "numeric"
normalized_counts <- sweep(count_matrix, 2L, size_factors$size_factor, "/")
normalized_output <- data.frame(peak_id = annotation$peak_id, normalized_counts, check.names = FALSE)
write_tsv(normalized_output, file.path(table_dir, "strict_candidate_normalized_counts_from_saved_size_factors.tsv"))

paired_donors <- sort(unique(metadata$donor[metadata$pairing_status == "complete_pair"]))
if (!identical(as.character(paired_donors), c("1001", "1003", "1004", "1008"))) {
  stop("Complete donor set differs from the protected four-pair design")
}

pair_summary <- lapply(seq_len(nrow(annotation)), function(i) {
  diffs <- vapply(paired_donors, function(donor_id) {
    mature_sample <- metadata$sample_id[metadata$donor == donor_id & metadata$cell_type == "Mature_NK"]
    memory_sample <- metadata$sample_id[metadata$donor == donor_id & metadata$cell_type == "Memory_NK"]
    if (length(mature_sample) != 1L || length(memory_sample) != 1L) {
      stop("Complete donor does not have exactly one Mature and one Memory sample: ", donor_id)
    }
    normalized_counts[i, memory_sample] - normalized_counts[i, mature_sample]
  }, FUN.VALUE = numeric(1))
  expected_sign <- sign(annotation$log2FoldChange_apeglm[i])
  consistent <- sign(diffs) == expected_sign
  data.frame(
    peak_id = annotation$peak_id[i],
    expected_apeglm_direction = annotation$apeglm_direction[i],
    paired_donors_consistent = sum(consistent),
    paired_donors_total = length(diffs),
    pair_consistency_fraction = mean(consistent),
    pair_differences_memory_minus_mature = paste(
      paste0(paired_donors, ":", formatC(diffs, format = "f", digits = 4)),
      collapse = ";"
    ),
    mean_normalized_mature_all_samples = mean(normalized_counts[i, metadata$cell_type == "Mature_NK"]),
    mean_normalized_memory_all_samples = mean(normalized_counts[i, metadata$cell_type == "Memory_NK"])
  )
})
pair_summary <- do.call(rbind, pair_summary)
write_tsv(pair_summary, file.path(table_dir, "strict_candidate_complete_pair_consistency.tsv"))

annotation <- merge(annotation, pair_summary, by = "peak_id", sort = FALSE)
annotation <- annotation[match(strict_annot_order$peak_id, annotation$peak_id), , drop = FALSE]
if (!identical(annotation$peak_id, strict_annot_order$peak_id)) {
  stop("Annotation order changed during pair-consistency merge")
}
write_tsv(annotation, file.path(table_dir, "strict_21_candidate_annotation_with_pair_consistency.tsv"))

gene_key <- ifelse(
  is.na(annotation$nearest_gene_symbol) | annotation$nearest_gene_symbol == "",
  annotation$nearest_gene_id_unversioned,
  annotation$nearest_gene_symbol
)
annotation$shortlist_gene_key <- gene_key
annotation$promoter_priority <- annotation$annotation_category == "Promoter"
annotation$absolute_apeglm_lfc <- abs(annotation$log2FoldChange_apeglm)
annotation$absolute_tss_distance <- abs(annotation$signed_distance_to_tss_bp)

rank_direction <- function(direction_label, n_take = 3L) {
  x <- annotation[annotation$apeglm_direction == direction_label, , drop = FALSE]
  x <- x[order(
    !x$promoter_priority,
    -x$paired_donors_consistent,
    -x$absolute_apeglm_lfc,
    x$rank_by_raw_pvalue,
    x$absolute_tss_distance,
    x$peak_id
  ), , drop = FALSE]
  x <- x[!duplicated(x$shortlist_gene_key), , drop = FALSE]
  head(x, n_take)
}

shortlist <- rbind(
  rank_direction("Memory_NK_more_accessible", 3L),
  rank_direction("Mature_NK_more_accessible", 3L)
)
if (nrow(shortlist) > 6L || anyDuplicated(shortlist$shortlist_gene_key)) {
  stop("Preliminary gene shortlist violates the six-gene or uniqueness rule")
}
shortlist$shortlist_status <- "quantitative_shortlist_pending_literature_review"
shortlist$selection_basis <- paste(
  "direction-balanced; promoter priority; four-pair consistency;",
  "absolute apeglm LFC; raw-p rank; TSS distance; one peak per nearest gene"
)
shortlist_columns <- c(
  "shortlist_gene_key", "nearest_gene_id", "nearest_gene_id_unversioned",
  "nearest_gene_symbol", "nearest_gene_type", "peak_id", "Chr", "Start", "End",
  "annotation_category", "signed_distance_to_tss_bp", "apeglm_direction",
  "log2FoldChange_apeglm", "log2FoldChange_unshrunken", "pvalue", "padj",
  "rank_by_raw_pvalue", "paired_donors_consistent", "paired_donors_total",
  "pair_consistency_fraction", "nearest_gene_ambiguous", "shortlist_status",
  "selection_basis"
)
write_tsv(shortlist[, shortlist_columns], file.path(table_dir, "preliminary_up_to_six_gene_shortlist.tsv"))

category_levels <- c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Downstream", "Distal/intergenic")
direction_labels <- c(
  Memory_NK_more_accessible = "Memory NK more accessible",
  Mature_NK_more_accessible = "Mature NK more accessible"
)
direction_colors <- c(
  "Memory NK more accessible" = "#2B6CB0",
  "Mature NK more accessible" = "#D97706"
)

category_counts <- as.data.frame(table(
  category = factor(annotation$annotation_category, levels = category_levels),
  direction = factor(direction_labels[annotation$apeglm_direction], levels = unname(direction_labels))
))
category_counts <- category_counts[category_counts$Freq > 0L, , drop = FALSE]
colnames(category_counts)[3L] <- "candidate_count"
category_totals <- aggregate(candidate_count ~ category, data = category_counts, sum)
write_tsv(category_counts, file.path(table_dir, "strict_candidate_genomic_category_counts_by_direction.tsv"))
write_tsv(category_totals, file.path(table_dir, "strict_candidate_genomic_category_totals.tsv"))

category_plot <- ggplot2::ggplot(
  category_counts,
  ggplot2::aes(x = category, y = candidate_count, fill = direction)
) +
  ggplot2::geom_col(width = 0.72) +
  ggplot2::geom_text(
    data = category_totals,
    ggplot2::aes(x = category, y = candidate_count, label = candidate_count),
    inherit.aes = FALSE, vjust = -0.35, fontface = "bold", size = 4.5
  ) +
  ggplot2::scale_fill_manual(values = direction_colors) +
  ggplot2::scale_y_continuous(
    limits = c(0, max(category_totals$candidate_count) + 1.4),
    breaks = seq(0, max(category_totals$candidate_count) + 1, by = 1),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::labs(
    title = "Genomic annotation of 21 strict primary nominal candidates",
    subtitle = "Promoter = +/-3 kb around a GENCODE v50 transcript TSS; zero FDR-significant DARs",
    x = NULL,
    y = "Candidate peaks",
    fill = "Direction"
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    panel.grid.major.x = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(angle = 25, hjust = 1),
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(figure_dir, "strict_candidate_genomic_annotation_categories.png"),
  category_plot, width = 10, height = 6.5, units = "in", dpi = 300, bg = "white"
)
ggplot2::ggsave(
  file.path(figure_dir, "strict_candidate_genomic_annotation_categories.pdf"),
  category_plot, width = 10, height = 6.5, units = "in", device = grDevices::cairo_pdf
)

distance_plot_data <- annotation
distance_plot_data$direction_label <- direction_labels[distance_plot_data$apeglm_direction]
distance_plot_data$signed_log10_distance <- sign(distance_plot_data$signed_distance_to_tss_bp) *
  log10(abs(distance_plot_data$signed_distance_to_tss_bp) + 1)
distance_plot_data$peak_label <- factor(
  distance_plot_data$peak_id,
  levels = distance_plot_data$peak_id[order(distance_plot_data$signed_distance_to_tss_bp)]
)
distance_ticks_bp <- c(-1e6, -1e5, -1e4, -3e3, 0, 3e3, 1e4, 1e5, 1e6)
distance_ticks_transformed <- sign(distance_ticks_bp) * log10(abs(distance_ticks_bp) + 1)
distance_tick_labels <- c("-1 Mb", "-100 kb", "-10 kb", "-3 kb", "0", "+3 kb", "+10 kb", "+100 kb", "+1 Mb")

distance_plot <- ggplot2::ggplot(
  distance_plot_data,
  ggplot2::aes(x = signed_log10_distance, y = peak_label, color = direction_label)
) +
  ggplot2::geom_vline(
    xintercept = sign(c(-3000, 0, 3000)) * log10(abs(c(-3000, 0, 3000)) + 1),
    linetype = c("dashed", "solid", "dashed"),
    color = c("#6B7280", "#111827", "#6B7280"), linewidth = c(0.6, 0.8, 0.6)
  ) +
  ggplot2::geom_point(size = 3.4, alpha = 0.9) +
  ggplot2::scale_color_manual(values = direction_colors) +
  ggplot2::scale_x_continuous(
    breaks = distance_ticks_transformed,
    labels = distance_tick_labels
  ) +
  ggplot2::labs(
    title = "Signed distance from each strict candidate to its nearest transcript TSS",
    subtitle = "Negative = upstream; positive = downstream; signed log scale; dashed guides mark +/-3 kb",
    x = "Signed distance to nearest TSS",
    y = "Candidate peak",
    color = "Direction"
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )

ggplot2::ggsave(
  file.path(figure_dir, "strict_candidate_signed_distance_to_nearest_tss.png"),
  distance_plot, width = 11, height = 8.2, units = "in", dpi = 300, bg = "white"
)
ggplot2::ggsave(
  file.path(figure_dir, "strict_candidate_signed_distance_to_nearest_tss.pdf"),
  distance_plot, width = 11, height = 8.2, units = "in", device = grDevices::cairo_pdf
)

annotation_audit <- data.frame(
  metric = c(
    "strict_candidates", "memory_more", "mature_more", "annotated_nearest_gene",
    "nearest_gene_ambiguous", "promoter_pm3kb", "unmapped_strict_candidates",
    "broader_homer_candidates", "primary_background_peaks",
    "background_contigs_without_gencode_features"
  ),
  value = c(
    nrow(annotation),
    sum(annotation$apeglm_direction == "Memory_NK_more_accessible"),
    sum(annotation$apeglm_direction == "Mature_NK_more_accessible"),
    sum(!is.na(annotation$nearest_gene_id) & annotation$nearest_gene_id != ""),
    sum(annotation$nearest_gene_ambiguous),
    sum(annotation$annotation_category == "Promoter"),
    0L,
    nrow(broader),
    nrow(ordinary_rows),
    length(setdiff(unique(to_gencode_seqname(ordinary_rows$Chr)), available_gencode_seqnames))
  )
)
write_tsv(annotation_audit, file.path(table_dir, "annotation_audit_summary.tsv"))

analysis_parameters <- data.frame(
  parameter = c(
    "execution_platform", "candidate_source", "strict_candidate_definition",
    "strict_candidate_count", "promoter_tss_window", "annotation_source",
    "annotation_assembly", "annotation_release", "annotation_level",
    "annotation_priority", "nearest_tss_distance_definition",
    "broader_homer_candidate_definition", "broader_homer_candidate_count",
    "motif_background", "motif_background_count", "planned_homer_window_bp",
    "deseq2_rerun", "apeglm_rerun", "maximum_threads", "interpretation_boundary"
  ),
  value = c(
    "Windows personal laptop",
    "protected primary 10-in-5 saved tables",
    "raw p-value <0.05 and absolute apeglm LFC >=1",
    "21",
    "-3000,+3000 bp around each transcript TSS",
    "GENCODE Release 50 comprehensive primary-assembly GTF",
    "GRCh38.p14",
    "Ensembl 116",
    "transcript",
    "Promoter > 5UTR > 3UTR > Exon > Intron > Downstream > Distal/intergenic",
    "signed strand-aware interval-boundary distance; zero when peak overlaps TSS",
    "raw p-value <0.05 and absolute ordinary LFC >=1",
    "1562",
    "all protected primary tested peaks",
    "59186",
    "200",
    "false",
    "false",
    "1 for this annotation attempt; project ceiling 2",
    "nominal exploratory hypotheses; zero FDR-significant DARs remains primary result"
  )
)
write_tsv(analysis_parameters, file.path(provenance_dir, "analysis_parameters.tsv"))

package_versions <- data.frame(
  package = c("R", required_packages),
  version = c(
    R.version.string,
    vapply(required_packages, function(pkg) as.character(utils::packageVersion(pkg)), character(1))
  )
)
write_tsv(package_versions, file.path(provenance_dir, "package_versions.tsv"))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) == 1L) sub("^--file=", "", script_arg) else NA_character_
execution_command <- paste(
  shQuote(file.path(R.home("bin"), "Rscript.exe")),
  "--vanilla",
  shQuote(normalizePath(script_path, winslash = "/", mustWork = TRUE)),
  shQuote(project_root),
  shQuote(gtf_path)
)
writeLines(execution_command, file.path(provenance_dir, "execution_command.txt"), useBytes = TRUE)

session_path <- file.path(provenance_dir, "session_info.txt")
session_connection <- file(session_path, open = "wt")
writeLines(capture.output(utils::sessionInfo()), session_connection)
close(session_connection)

metrics <- data.frame(
  metric = c(
    "strict_candidates", "memory_more", "mature_more", "broader_homer_candidates",
    "background_peaks", "annotation_rows", "shortlisted_genes",
    "all_four_pairs_consistent_candidates", "promoter_candidates"
  ),
  value = c(
    21L, 11L, 10L, 1562L, 59186L, nrow(annotation), nrow(shortlist),
    sum(pair_summary$paired_donors_consistent == 4L),
    sum(annotation$annotation_category == "Promoter")
  )
)
write_tsv(metrics, file.path(provenance_dir, "analysis_metrics.tsv"))

status <- data.frame(
  status_key = c(
    "WEEK3_ANNOTATION_ATTEMPT2_EXECUTION",
    "WEEK3_ANNOTATION_ATTEMPT2_INPUT_HASHES",
    "WEEK3_ANNOTATION_ATTEMPT2_CANDIDATE_COUNTS",
    "WEEK3_ANNOTATION_ATTEMPT2_BED_COORDINATES",
    "WEEK3_ANNOTATION_ATTEMPT2_ANNOTATION_ROWS",
    "WEEK3_ANNOTATION_ATTEMPT2_INTERPRETATION"
  ),
  status_value = c(
    "COMPLETED",
    "PASS_ALL_7",
    "PASS_21_STRICT_1562_BROADER_59186_BACKGROUND",
    "PASS_ZERO_BASED_HALF_OPEN_ROUNDTRIP",
    "PASS_21_OF_21_NEAREST_TSS_ASSIGNED",
    "EXPLORATORY_NOMINAL_CANDIDATES_ZERO_FDR_DARS_PRESERVED"
  )
)
write_tsv(status, file.path(provenance_dir, "analysis_status.tsv"))

readme_lines <- c(
  "# Week 3 downstream annotation - attempt 2",
  "",
  "This Windows-local attempt reuses protected Week 2 tables and does not rerun DESeq2 or apeglm.",
  "",
  "- Strict set: 21 primary nominal candidates (raw p <0.05 and absolute apeglm LFC >=1).",
  "- Annotation: GENCODE Release 50 comprehensive primary assembly, GRCh38.p14 / Ensembl 116.",
  "- Promoter: plus or minus 3 kb around a transcript TSS.",
  "- Broader HOMER sensitivity BED: 1,562 primary raw-p plus ordinary-LFC candidates.",
  "- Custom motif background BED: all 59,186 primary tested peaks.",
  "",
  "All peak-gene links and the preliminary six-gene shortlist are exploratory hypotheses.",
  "The protected primary result remains zero FDR-significant DARs."
)
writeLines(readme_lines, file.path(result_root, "README.md"), useBytes = TRUE)

manifest_files <- list.files(result_root, recursive = TRUE, full.names = TRUE)
manifest_files <- manifest_files[file.info(manifest_files)$isdir %in% FALSE]
manifest <- data.frame(
  relative_path = substring(normalizePath(manifest_files, winslash = "/"), nchar(normalizePath(result_root, winslash = "/")) + 2L),
  bytes = file.info(manifest_files)$size,
  sha256 = vapply(manifest_files, sha256_file, FUN.VALUE = character(1))
)
manifest <- manifest[order(manifest$relative_path), , drop = FALSE]
write_tsv(manifest, file.path(provenance_dir, "output_sha256_manifest.tsv"))

message("Week 3 annotation attempt 2 completed successfully")
