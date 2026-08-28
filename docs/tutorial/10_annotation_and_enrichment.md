# 10. Annotation and enrichment

## Annotation

Use an annotation release compatible with the reference assembly and chromosome
names. Promoter, intronic, exonic, and distal categories depend on declared
rules. A nearest-gene or nearest-TSS label is a proximity hypothesis, not proof
that a peak regulates that gene.

## Enrichment background

The background should represent regions or genes that could have been selected
by the analysis. Using the whole genome when only a filtered subset was tested
can bias enrichment. Record peak-to-gene mapping, duplicate-gene handling,
database version, correction method, and the exact target and background sizes.

## Motif analysis

Use an explicit target set and matched accessible-region background. Motif
enrichment indicates sequence-pattern enrichment; it does not establish that a
transcription factor is expressed, bound, active, or causal.

## When there are no FDR-significant regions

Lead with the zero-FDR primary result. Optional raw-p-value or effect-size
candidate sets must be labelled nominal and exploratory everywhere—in filenames,
figures, captions, tables, and prose. They must not replace the primary result.

## Downstream gate

- target definition is frozen and reported;
- background derives from the tested universe;
- genome and annotation releases match;
- multiple testing is handled within enrichment tools;
- proximity, enrichment, and motif findings use non-causal language.
