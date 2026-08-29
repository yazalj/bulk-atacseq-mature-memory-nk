# 02. Samples, metadata, and statistical design

## Build the sample registry first

Create one row per biological library using
[`samples.template.tsv`](../../config/templates/samples.template.tsv). Use
stable sample identifiers and explicit columns for donor or subject, condition,
batch, treatment, and read mates. Do not infer biology solely from filenames.

For paired-end data, verify that every sample has exactly one R1 and one R2,
both belong to the same library, and read counts are plausible. Record file
checksums in the input manifest.

## Copy-and-adapt metadata checks

Fill `config/local/samples.tsv` first. The final BAM column names the output
that will be created in Chapter 6 and counted in Chapter 8. From the project
root, check the table and the input files:

```bash
set -euo pipefail

SAMPLES="config/local/samples.tsv"
test -s "${SAMPLES}"

awk -F '\t' '
  NR == 1 {
    for (i = 1; i <= NF; i++) header[$i] = i
    required[1] = "sample_id"
    required[2] = "donor"
    required[3] = "condition"
    required[4] = "read1_fastq"
    required[5] = "read2_fastq"
    required[6] = "final_bam"
    for (i = 1; i <= 6; i++) {
      if (!(required[i] in header)) {
        print "Missing column: " required[i] > "/dev/stderr"
        exit 1
      }
    }
    next
  }
  {
    if ($header["sample_id"] == "") exit 1
    if (seen[$header["sample_id"]]++) {
      print "Duplicate sample_id: " $header["sample_id"] > "/dev/stderr"
      exit 1
    }
  }
  END { if (NR < 2) exit 1 }
' "${SAMPLES}"

awk -F '\t' '
  NR == 1 { for (i = 1; i <= NF; i++) header[$i] = i; next }
  { print $header["sample_id"] "\t" $header["read1_fastq"] "\t" $header["read2_fastq"] }
' "${SAMPLES}" |
while IFS=$'\t' read -r sample_id read1 read2; do
  test -s "${read1}" || { echo "Missing R1: ${read1}" >&2; exit 1; }
  test -s "${read2}" || { echo "Missing R2: ${read2}" >&2; exit 1; }
  gzip -t "${read1}"
  gzip -t "${read2}"
done < "${SAMPLES}"

echo "Metadata and FASTQ presence checks: PASS"
```

Create checksums for the immutable inputs and paste them into the input
manifest. Repeat the command for every read and reference file:

```bash
sha256sum /path/to/sample_01_R1.fastq.gz /path/to/sample_01_R2.fastq.gz
sha256sum /path/to/reference.fa /path/to/annotation.gtf /path/to/blacklist.bed
```

## Define the comparison

Write the intended statistical formula before counting. Examples include:

```r
~ condition
~ batch + condition
~ donor + condition
```

The correct design depends on the experiment. A paired-donor term is valuable
only when donor identities and condition assignments support it. Confounded
designs cannot be repaired by adding more terms. Build the model matrix and
confirm that it is full rank before fitting DESeq2.

After replacing the example design with the intended formula, this R block
performs an early rank check:

```r
metadata <- read.delim("config/local/samples.tsv", check.names = FALSE)
metadata$donor <- factor(metadata$donor)
metadata$condition <- factor(metadata$condition)

design_formula <- ~ donor + condition  # CHANGE_ME if the study needs another design
model_matrix <- model.matrix(design_formula, data = metadata)

if (qr(model_matrix)$rank != ncol(model_matrix)) {
  stop("The proposed design matrix is not full rank. Resolve confounding first.")
}

print(levels(metadata$condition))
print(colnames(model_matrix))
```

## Metadata gate

- sample IDs are unique;
- R1/R2 mates are complete and correctly paired;
- count-matrix column order will be explicit;
- factor levels and contrast direction are written in plain language;
- biological replicates, not sequencing lanes, define replication;
- batch and condition are not perfectly confounded;
- exclusions are justified before inspecting the final contrast.
