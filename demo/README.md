# Synthetic count-matrix smoke test

This tiny, fully synthetic example checks that the repository can validate an
ordered paired-donor sample table, read a nonnegative integer count matrix,
apply a declared low-count rule, and reproduce known output files. It uses only
Python's standard library and normally finishes in less than a second.

It **does not** run read QC, alignment, BAM processing, peak calling,
featureCounts, DESeq2, or biological interpretation. Its descriptive means and
log2 fold changes are not p-values, adjusted p-values, or differential-
accessibility results.

## Run the read-only check

From the repository root, paste:

```bash
python3 scripts/run_synthetic_demo.py --check
```

Expected terminal message:

```text
Synthetic demo: PASS (5 of 6 peaks retained)
```

The command uses a temporary directory and compares all generated files with
the tracked files under [`expected/`](expected/).

## Inspect a new output attempt

To retain the generated files locally, choose a new attempt directory:

```bash
python3 scripts/run_synthetic_demo.py --output-dir demo/run_attempt1
```

The script refuses to reuse an existing output directory. The generated
`retained_counts.tsv`, `descriptive_summary.tsv`, and `run_summary.json` can be
compared with [`expected/`](expected/). Local `demo/run_attempt*/` directories
are ignored by Git.

## Deliberate validation behavior

Try the test only after making a separate copy of the inputs if you want to
experiment. It stops when sample IDs are duplicated or reordered, a donor lacks
one condition, counts are negative or non-integer, peak IDs repeat, or the
filter retains no peaks. This behavior illustrates the validation style used
throughout the full tutorial and the AI-agent contract.
