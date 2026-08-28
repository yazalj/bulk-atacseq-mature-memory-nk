---
name: bulk-atacseq-agent
description: Plan, validate, or advance a human-supervised bulk paired-end ATAC-seq analysis using this repository. Use for AI-assisted dataset readiness checks, configuration, staged execution, validation, or reporting. Do not use for single-cell ATAC-seq, multiome data, or autonomous end-to-end execution.
---

# Bulk ATAC-seq Agent

Use the repository's provider-neutral contract as the source of truth.

## Begin every task

1. Read `../../../AGENTS.md`.
2. Read `../../../agent/README.md` and
   `../../../agent/workflow_contract.yaml`.
3. Determine whether the request is `case_study_audit`, `new_study_dry_run`, or
   `new_study_execution`. Default a new dataset to dry-run mode.
4. Read only the tutorial chapters relevant to the current stage.
5. Inspect the supplied sample, parameter, manifest, and status files.
6. Run `../../../scripts/validate_agent_configuration.py` before proposing
   biological computation.

## Preserve human control

- Ask for a missing scientific decision; never synthesize one from filenames or
  convenience.
- Treat approval as stage-specific. Do not interpret a broad request as
  permission to bypass failed gates or unresolved decisions.
- In execution mode, run only the one approved stage, validate it, update the
  project-specific status record, and stop before the next stage.
- Keep raw and reference inputs immutable and use non-overwriting attempt paths.
- Do not install software, change environments, upload data, publish results,
  or use credentials without explicit authorization.

## Use the implementation honestly

The files under `../../../workflow/` are case-study implementations with locked
NK-cell sample and result assertions. Do not run them unchanged on a new study.
Adapt assumptions and validation together after the new design is approved.

Use adjusted-p-value-controlled results for primary claims. Keep nominal
candidates, nearest-gene assignments, enrichment, motifs, and locus views
explicitly exploratory when the primary evidence does not support discovery.

## Finish a turn

Report the active mode, configuration-validation result, completed evidence,
unresolved decisions, current status, and the next single approvable stage. Do
not claim completion without the contract's validation evidence.
