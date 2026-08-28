# 08. Consensus peaks and fragment counting

## Define the feature universe

A differential analysis needs the same genomic features for every sample.
Declare how replicate peak sets become a consensus universe. Options include a
simple union, a reproducibility-supported set, or a study-specific consensus
rule. The choice affects both sensitivity and the multiple-testing burden.

The internship merged peaks within each condition, combined the condition
unions, and merged overlapping or book-ended intervals with `bedtools merge
-d 0`. Singleton-supported regions were retained. This created 112,759 regions.

The worked implementation is
[`workflow/week2_consensus_counts/`](../../workflow/week2_consensus_counts/).
Its safety checks are intentionally locked to the internship's ten samples, so
a new study must revise the manifest and expected dimensions deliberately.

## Count paired fragments

Convert the BED universe to the annotation format required by the counter while
preserving coordinate semantics. The internship used one-based inclusive SAF
coordinates and featureCounts:

```text
-F SAF -p --countReadPairs -B -C -s 0 -T 2
```

## Counting gate

- consensus intervals are valid, sorted, and non-overlapping;
- coordinate conversion preserves interval lengths exactly;
- BAM order is explicit and matches the planned metadata order;
- counts are nonnegative integers and represent fragments, not individual mates;
- column sums agree with the counter's assignment summary;
- the clean matrix contains no silent filtering, normalization, or reordering.
