# 11. Interpretation and reporting

## Report in this order

1. Experimental design, samples, exclusions, and reference build.
2. Processing and QC decisions with versions and parameters.
3. Consensus definition, count unit, filter, model, reference level, and
   contrast.
4. Number of regions entering and passing each gate.
5. Primary multiple-testing-controlled result and effect-size threshold.
6. Diagnostic figures and limitations.
7. Clearly separated exploratory annotation or enrichment.

## Language for common outcomes

Prefer “no region met the prespecified threshold” over “there was no biological
difference.” A non-significant result can reflect limited power, heterogeneity,
measurement noise, small effects, or a genuinely weak contrast.

Do not describe nearest genes as targets, motifs as active transcription
factors, visual coverage differences as normalized group effects, or nominal
p-values as discoveries.

## Reproducibility package

Preserve metadata, input checksums, parameters, scripts, software versions,
commands, validation summaries, compact result tables, figure sources, and a
limitations document. Large or controlled data should remain in an appropriate
archive rather than Git. Make public-data accessions and regeneration steps
discoverable.

## Copy-and-adapt provenance commands

Run these commands from the project root after activating the same environment
used for analysis. The `environment/local/` directory is ignored by this
repository so machine-specific package URLs are not published accidentally.

```bash
set -euo pipefail

mkdir -p environment/local logs/reporting
conda list --explicit > environment/local/conda_explicit_spec.txt
conda env export > environment/local/conda_environment_resolved.yml
Rscript -e 'sessionInfo()' > environment/local/R_sessionInfo.txt

sha256sum \
  config/local/samples.tsv \
  config/local/input_manifest.tsv \
  config/local/analysis_parameters.yaml \
  results/consensus/consensus.bed \
  results/counts/consensus_counts.clean.tsv \
  results/differential_accessibility/all_tested_regions.tsv \
  > logs/reporting/key_output_sha256.txt

date '+%Y-%m-%dT%H:%M:%S%z' > logs/reporting/reporting_timestamp.txt
```

If a listed output does not exist because that optional stage was not run,
remove only that path from the checksum command and state why it is absent.
Keep the exact commands used for each sample and stage alongside these records.

The worked example's public reporting choices can be reviewed in
[`docs/METHODS.md`](../METHODS.md), [`docs/LIMITATIONS.md`](../LIMITATIONS.md),
and [`docs/REPRODUCIBILITY.md`](../REPRODUCIBILITY.md).
