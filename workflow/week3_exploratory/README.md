# Exploratory Week 3 scripts

[← Previous: Differential accessibility](../week2_differential_accessibility/README.md) · [Workflow contents](../README.md)

These scripts operate on nominal candidate sets after the primary analysis
returned zero FDR-significant DARs.

- `07_annotate_nominal_candidates.R` annotates the strict 21-region nominal set
  with GENCODE Release 50 and nearest-transcript-TSS hypotheses.
- `08_prepare_gprofiler_inputs.py` prepares direction-aware gene queries and a
  custom background for manual g:Profiler submission.
- `09_run_homer_nominal_sensitivity.sh` records the one-worker HOMER sensitivity
  analysis using 1,562 nominal target regions, the 59,186-region primary
  background universe, and central 200-bp windows.

The annotation and g:Profiler preparation scripts require the full omitted
peak-level result tables plus verified GENCODE resources. The HOMER script
requires explicit `PROJECT_ROOT`, `HOMER_ROOT`, `GENOME_FASTA`, and
`GENOME_FASTA_GZ` environment variables. Its paths must use the syntax expected
by the shell in which HOMER is installed.

None of these analyses establishes a significant DAR, target gene,
transcription-factor activity, or regulatory mechanism.

---

[← Previous: Differential accessibility](../week2_differential_accessibility/README.md) · [Workflow contents](../README.md)
