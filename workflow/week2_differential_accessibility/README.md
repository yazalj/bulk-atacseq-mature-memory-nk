# Differential-accessibility scripts

`04_run_primary_deseq2.R` reproduces the fixed primary analysis from a clean
count matrix, peak annotation, and ordered metadata. It uses the 10-in-5 filter,
design `~ donor + cell_type`, Mature NK reference, Memory-relative-to-Mature
contrast, adjusted-p-value threshold 0.05, absolute-LFC threshold 1, and apeglm
shrinkage.

The script can install missing R packages into the active user library. Review
that behavior and the recorded versions before running it in a new environment.
It writes to a new `public_reproduction_run` directory and refuses to overwrite
an existing run.

`05_summarize_nominal_sensitivity.R` and
`06_plot_nominal_sensitivity_volcano.R` require the full saved result tables for
the primary 10-in-5, 10-in-3, and no-explicit-prefilter analyses. Those large
tables are not included in this source-oriented repository. The scripts summarize
existing tables; highlighted raw-p-value regions are nominal, not FDR-significant.
