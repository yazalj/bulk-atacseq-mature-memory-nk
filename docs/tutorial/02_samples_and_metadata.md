# 02. Samples, metadata, and statistical design

## Build the sample registry first

Create one row per biological library using
[`samples.template.tsv`](../../config/templates/samples.template.tsv). Use
stable sample identifiers and explicit columns for donor or subject, condition,
batch, treatment, and read mates. Do not infer biology solely from filenames.

For paired-end data, verify that every sample has exactly one R1 and one R2,
both belong to the same library, and read counts are plausible. Record file
checksums in the input manifest.

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

## Metadata gate

- sample IDs are unique;
- R1/R2 mates are complete and correctly paired;
- count-matrix column order will be explicit;
- factor levels and contrast direction are written in plain language;
- biological replicates, not sequencing lanes, define replication;
- batch and condition are not perfectly confounded;
- exclusions are justified before inspecting the final contrast.
