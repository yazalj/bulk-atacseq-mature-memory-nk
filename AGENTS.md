# Repository instructions for AI agents

## Scope

This repository contains two related resources:

1. a preserved public case study of a completed Mature-versus-Memory NK-cell
   bulk ATAC-seq internship analysis; and
2. a reusable, human-supervised guide for new bulk, paired-end ATAC-seq studies.

Do not treat this repository as a universal or fully autonomous pipeline. It
does not cover single-cell ATAC-seq, multiome assays, clinical-grade processing,
or every organism and library design.

## Required reading

Before proposing analysis commands or changing analysis state, read:

1. `agent/README.md`;
2. `agent/workflow_contract.yaml`;
3. `docs/tutorial/README.md`;
4. `docs/tutorial/02_samples_and_metadata.md`;
5. the tutorial chapter for the proposed stage; and
6. the project-specific sample, parameter, and status files supplied by the
   user.

For the completed case study, also read `docs/REPRODUCIBILITY.md` and
`config/analysis_parameters.yaml`.

## Operating modes

Use exactly one mode and state it in the analysis-status file:

- `case_study_audit`: inspect the published example without rerunning it.
- `new_study_dry_run`: validate inputs and prepare a proposed staged plan; do
  not run biological computation.
- `new_study_execution`: execute only the single stage explicitly authorized by
  the user after its dependencies and human decisions are complete.

Default to `new_study_dry_run` for a new dataset unless the user clearly
authorizes execution.

## Non-negotiable safeguards

- Treat raw reads, reference resources, and supplied source data as immutable.
- Never infer sample identity, mates, organism, genome build, adapter sequence,
  blacklist, reference level, contrast, batch structure, or biological pairing
  from filenames alone.
- Never invent a scientific parameter to make a command runnable. Record the
  unresolved decision and stop for human input.
- Never overwrite an input, completed stage, or prior attempt. Use a new,
  clearly named output directory.
- Do not install or modify software environments, upload data, publish results,
  or use credentials without explicit authorization.
- Do not run computationally intensive stages concurrently. Default to no more
  than two total CPU threads until the user records a different authorized
  limit.
- Preserve paired-end relationships. Blacklist removal must be fragment-aware:
  if either mate is excluded, both alignments of that fragment must be removed.
- Validate sample order before counting and again before statistical modelling.
- Use adjusted p-values for primary genome-wide claims. Never silently promote
  nominal raw-p-value candidates to discoveries.
- Treat nearest-gene, enrichment, motif, and locus-view results as hypotheses
  unless supported by the declared primary evidence and independent validation.
- Do not claim completion from command exit status alone; run the required
  stage validations and record their evidence.

## Human decisions and authorization

The human owner must approve the study scope, reference resources, read-
processing policy, alignment/filtering policy, peak-calling policy, consensus
strategy, statistical design and contrast, and downstream-analysis scope.

Approval for one stage does not authorize later stages. When a required
decision is missing, update the status record to `awaiting_human` and ask a
focused question. Do not continue by assumption.

## Configuration and state

- Start from the files under `config/templates/` and `agent/`.
- Prefer the JSON parameter format for deterministic agent validation. The YAML
  template remains available for human-facing configuration.
- Run `scripts/validate_agent_configuration.py` before proposing execution.
- Keep runtime state in a project-specific copy of
  `agent/analysis_status.template.json`; do not edit the template itself.
- Record input identity, command, tool path and version, outputs, checks, warnings,
  and approval evidence for each completed stage.

## Case-study scripts

The scripts under `workflow/` contain deliberate checks for the recorded
internship design, including sample identities, dimensions, factor levels, and
thresholds. They are reference implementations, not drop-in scripts for an
arbitrary dataset. Adapt their assumptions only after the new study design is
approved, and update validation assertions together with the implementation.

Do not rerun or overwrite the published case-study results merely because the
scripts are present.

## Repository maintenance

When changing repository-facing files:

- preserve the zero-FDR-first interpretation of the case study;
- keep relative Markdown links valid;
- keep credentials, private paths, controlled data, and large sequencing files
  outside Git;
- update `MANIFEST.sha256` mechanically after the final content changes;
- run relevant validators and syntax checks; and
- do not commit or push unless the user requested publication.

## Completion criteria

A stage is complete only when its declared outputs exist, required checks pass,
the status record contains the evidence, and no stop condition remains. A full
analysis is complete only when reporting and limitations are also recorded.
