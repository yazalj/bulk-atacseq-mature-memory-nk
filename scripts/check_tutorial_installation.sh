#!/usr/bin/env bash

set -uo pipefail

failures=0

check_program() {
  local program="$1"
  if command -v "${program}" >/dev/null 2>&1; then
    printf 'PASS  %-24s %s\n' "${program}" "$(command -v "${program}")"
  else
    printf 'MISS  %-24s required command not found\n' "${program}"
    failures=$((failures + 1))
  fi
}

printf '%s\n' 'Bulk ATAC-seq tutorial installation check'
printf '%s\n' '========================================='

for program in fastqc multiqc cutadapt bowtie2 bowtie2-build samtools \
  picard bedtools macs3 featureCounts R Rscript; do
  check_program "${program}"
done

for program in annotatePeaks.pl findMotifsGenome.pl; do
  if command -v "${program}" >/dev/null 2>&1; then
    printf 'PASS  %-24s %s\n' "${program}" "$(command -v "${program}")"
  else
    printf 'INFO  %-24s optional HOMER command not found\n' "${program}"
  fi
done

if command -v Rscript >/dev/null 2>&1; then
  if Rscript -e '
required <- c("DESeq2", "apeglm", "vsn", "ggplot2", "pheatmap", "data.table", "gprofiler2")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  message("Missing R packages: ", paste(missing, collapse = ", "))
  quit(status = 1)
}
cat("PASS  Required R packages are available.\n")
'; then
    :
  else
    failures=$((failures + 1))
  fi
fi

if (( failures > 0 )); then
  printf '\nInstallation verification: FAIL (%d required check(s) missing)\n' "${failures}" >&2
  exit 1
fi

printf '\n%s\n' 'Installation verification: PASS'
