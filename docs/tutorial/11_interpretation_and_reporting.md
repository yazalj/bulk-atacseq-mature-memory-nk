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

The worked example's public reporting choices can be reviewed in
[`docs/METHODS.md`](../METHODS.md), [`docs/LIMITATIONS.md`](../LIMITATIONS.md),
and [`docs/REPRODUCIBILITY.md`](../REPRODUCIBILITY.md).
