# Consensus peaks and paired-fragment counts

These three scripts are concise public-release adaptations. They preserve the
sample order, merge rules, fragment-counting options, validation checks, and
overwrite protection of the validated workflow, but they were not executed to
create a new public result set.

Required environment variables:

```bash
export PROJECT_ROOT=/path/to/bulk-atacseq-mature-memory-nk
export INPUT_ROOT=/path/to/validated_bam_bai_and_broadpeak_inputs
export TOOL_DIR=/optional/path/to/analysis/bin
export EXPECTED_HOST=/optional/short_hostname
```

Run from `PROJECT_ROOT`, one stage at a time:

```bash
ATTEMPT=public_reproduction bash workflow/week2_consensus_counts/01_build_consensus_peaks.sh
ATTEMPT=public_reproduction bash workflow/week2_consensus_counts/02_count_fragments.sh
ATTEMPT=public_reproduction python3 workflow/week2_consensus_counts/03_create_clean_count_matrix.py
```

The expected input basenames and fixed biological order are defined in
`config/week2_samples.tsv`. The scripts refuse to overwrite an existing attempt.
