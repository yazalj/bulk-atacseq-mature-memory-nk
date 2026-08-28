# Workflow scripts

These scripts are organized by analytical stage. They are public-release
adaptations of preserved project scripts, not newly executed analyses.

- `week1_training/`: method record for the separate chromosome-22 exercise.
- `week2_consensus_counts/`: consensus construction, paired-fragment counting,
  and clean-matrix creation.
- `week2_differential_accessibility/`: primary donor-aware DESeq2 analysis and
  nominal sensitivity summaries.
- `week3_exploratory/`: nominal-candidate annotation, g:Profiler input
  preparation, and HOMER sensitivity analysis.

The scripts intentionally refuse to overwrite existing output locations. Review
all paths and software requirements before execution.
