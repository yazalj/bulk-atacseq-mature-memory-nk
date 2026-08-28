# Limitations and interpretation boundaries

- The primary donor-aware analysis detected zero regions meeting the declared
  FDR and effect-size thresholds. This is the controlling result.
- The dataset contains four complete donor pairs and two Memory-only donors.
  The model adjusts for donor, but the design remains small and unbalanced.
- Peak-universe construction used condition-level unions and retained peaks
  supported by a single sample. Conclusions depend on this declared universe.
- Count-assignment rates varied substantially between samples and are technical
  context, not a basis for post hoc exclusion.
- Raw-p-value candidate sets are nominal. They do not control the false
  discovery rate.
- Apeglm shrinkage reduced many ordinary effects, emphasizing uncertainty in
  low-information estimates.
- Nearest-gene labels do not demonstrate enhancer-gene regulation or gene-
  expression changes.
- g:Profiler terms depend on nominal region selection and nearest-gene mapping.
  They do not prove pathway activation.
- A HOMER motif label describes sequence similarity and enrichment relative to
  a background. It does not demonstrate transcription-factor expression,
  binding, occupancy, or causality.
- IGV tracks were independently autoscaled and were used for qualitative visual
  QC, not normalized group comparison.
- Week 1 used one Immature NK sample and a chromosome-22 teaching subset. Its QC
  cannot be generalized to the full ten-sample comparison.
- No footprinting analysis was completed; the repository makes no footprinting
  or differential-binding claim.
