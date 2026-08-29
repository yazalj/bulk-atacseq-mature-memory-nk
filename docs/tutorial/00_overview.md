# 00. Overview and decision points

[Tutorial contents](README.md) · [Next: Resource planning →](RESOURCE_PLANNING.md)

## Objective

Convert paired-end bulk ATAC-seq reads into a quality-controlled region-by-
sample fragment-count matrix, test a prespecified accessibility contrast, and
report the result with appropriate multiple-testing control and limitations.

## Decisions that must not be copied blindly

Before running commands, decide and record:

- organism, genome assembly, chromosome naming convention, and reference
  release;
- library layout, adapter sequences, and whether mitochondrial or alternate
  contigs are retained;
- biological unit, replicate structure, batches, paired donors, covariates,
  reference level, and contrast;
- mapping-quality, duplicate, blacklist, peak-calling, and low-count rules;
- whether peaks are narrow or broad and how the consensus universe is built;
- primary adjusted-p-value and effect-size thresholds;
- the universe used for enrichment and any explicitly exploratory analyses.

The internship used GRCh38, paired-end libraries, a donor-aware model, and a
10-count-in-5-samples filter because those choices matched that dataset and the
agreed analysis plan. They are an example, not universal defaults.

## Reproducibility rule

For every stage, preserve four things: immutable input identity, the exact
command and software version, validation evidence, and a new non-overwriting
output location. Separate observed results from interpretation.

---

[Tutorial contents](README.md) · [Next: Resource planning →](RESOURCE_PLANNING.md)
