# Reproducibility and provenance

## Status of this repository

This repository was assembled after completion of the internship. It is a
curated release snapshot, not the original Git execution history. The protected
analysis outputs were not rerun while preparing it.

The public scripts are adapted copies of validated project scripts. Adaptations
remove personal paths, rename tutor-facing output locations, and expose local
paths through arguments or environment variables. Scientific definitions,
sample order, model design, thresholds, and thread limits are retained.

## Suggested execution order

1. Review and adapt `config/input_manifest.example.tsv` without committing
   private absolute paths.
2. Verify sample identities and order against `config/samples_public.tsv`.
3. Export `PROJECT_ROOT`, `INPUT_ROOT`, and optionally `TOOL_DIR` and
   `EXPECTED_HOST`.
4. Run the three consensus/counting scripts sequentially.
5. Run `04_run_primary_deseq2.R` from the repository root after the count
   matrix, annotation, and ordered metadata are present.
6. Run sensitivity and Week 3 scripts only if their required full result tables
   and reference resources have been regenerated.

No analytical step should be run concurrently. The recorded ceiling was two
total CPU threads and one intensive process at a time.

## Important scope distinction

The small tables in `results/summary_tables/` support review of the reported
summary. They are not sufficient to rerun DESeq2, annotation, HOMER, or
g:Profiler. Full count and peak-level result tables are omitted because this is
a source-oriented public repository rather than a data archive.

## Included evidence

- public sample registry;
- declared analysis parameters;
- public-release workflow scripts;
- small validated numerical summaries;
- selected validated figures;
- recorded software versions;
- an inclusion checksum manifest generated for this release.

## Excluded evidence

- raw and processed sequencing data;
- reference files and software installations;
- large matrices and peak-level result tables;
- raw runtime logs containing workstation or server paths;
- internal execution handoffs and tutor instructions;
- failed-attempt, rendering, and build histories;
- university-private or course-supplied materials.
