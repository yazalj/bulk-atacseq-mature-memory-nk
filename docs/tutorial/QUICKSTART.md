# Quick-start checklist

[Tutorial contents](README.md) · [Start the full tutorial: Overview →](00_overview.md)

Use this page as a navigation aid, not as a substitute for the chapter-level
checks.

1. Follow the [installation and setup chapter](01_requirements_and_project_setup.md),
   activate the tutorial environment, and verify every required executable and
   R package before downloading or processing data.
2. State the biological comparison, experimental unit, replicate structure,
   likely nuisance variables, organism, and reference build.
3. Copy the files described in the [`config/templates/` guide](../../config/templates/README.md) and replace
   every `REVIEW_ME`, placeholder path, and example sample.
4. Record checksums for reads and reference resources. Confirm that FASTA,
   indexes, annotation, chromosome names, and blacklist all use the same build.
5. Run FastQC on raw mates and summarize only the intended reports with
   MultiQC. Decide whether trimming is justified by the evidence.
6. If trimming, process both mates together and rerun read-level QC.
7. Align paired reads to the verified reference. Record the executable and
   version, command, exit status, and mapping summary.
8. Apply declared mapping-quality, flag, duplicate, unwanted-contig, and
   pair-preserving blacklist rules. Validate BAM integrity, pairing, fragment-
   length periodicity, mitochondrial fraction, and TSS enrichment where the
   required annotation is suitable.
9. Call peaks independently for every biological sample using paired-fragment
   mode. Check peak counts, peak widths, and FRiP alongside the earlier QC
   metrics when reviewing outliers.
10. Construct a declared consensus universe and count paired fragments in an
   explicitly recorded sample order.
11. Prove that count-matrix columns exactly match the metadata rows. Check that
    the proposed model matrix is full rank.
12. Choose and record the low-count rule before testing. Fit the model with the
    intended reference level and contrast; apply multiple-testing correction.
13. Inspect size factors, PCA, sample distances, dispersion behavior, and
    effect estimates. Treat these as diagnostics, not significance tests.
14. Annotate and enrich only with a suitable tested-region background. If the
    primary analysis has no FDR-significant regions, report that result first.
15. Preserve parameters, versions, logs, checksums, small summaries, figures,
    and limitations. Do not publish controlled data or credentials.

Proceed to the next stage only after the current stage's expected outputs and
validation gate have passed.

---

[Tutorial contents](README.md) · [Start the full tutorial: Overview →](00_overview.md)
