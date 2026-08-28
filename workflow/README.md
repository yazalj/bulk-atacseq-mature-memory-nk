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

## Relationship to the reusable tutorial

The [step-by-step tutorial](../docs/tutorial/README.md) explains the decisions,
inputs, outputs, and validation gates around each stage. The scripts in this
directory remain the worked internship implementation and contain deliberate
case-study checks, including the recorded sample order and design. They are not
advertised as a universal push-button pipeline.

| Tutorial stage | Case-study implementation |
|---|---|
| Raw-read QC through peak calling | [`week1_training/README.md`](week1_training/README.md) and [`docs/METHODS.md`](../docs/METHODS.md) |
| Consensus regions and fragment counting | [`week2_consensus_counts/`](week2_consensus_counts/) |
| Differential accessibility | [`week2_differential_accessibility/`](week2_differential_accessibility/) |
| Annotation and enrichment | [`week3_exploratory/`](week3_exploratory/) |

For a new dataset, begin with the blank files under
[`config/templates/`](../config/templates/) and revise the scripts only after
the new metadata, reference resources, thresholds, and statistical design have
been justified.
