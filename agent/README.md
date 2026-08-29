# AI-agent interface

This directory provides a provider-neutral contract for using an AI agent to
prepare, validate, and—only with staged human approval—execute a bulk,
paired-end ATAC-seq analysis based on this repository's tutorial.

An AI agent is an orchestration and review layer. Deterministic command-line
tools perform the bioinformatics, and deterministic validators check important
invariants. The agent must not replace scientific decisions with guesses.

## What is included

- [`workflow_contract.yaml`](workflow_contract.yaml): ordered stages,
  dependencies, human decisions, expected evidence, and stop conditions.
- [`analysis_status.template.json`](analysis_status.template.json): durable
  machine-readable state for one analysis.
- JSON Schemas for the [sample registry](schemas/samples.schema.json),
  [parameter document](schemas/parameters.schema.json), and
  [analysis status](schemas/status.schema.json).
- [`examples/`](examples/README.md): a small configuration-only example with no
  sequencing data.
- [`example_request.md`](example_request.md): a reusable prompt for beginning a
  safe dry run.
- [`../scripts/validate_agent_configuration.py`](../scripts/validate_agent_configuration.py):
  standard-library validation for TSV, JSON, and status consistency.

The root [`AGENTS.md`](../AGENTS.md) supplies repository instructions to agents
that support that convention. The repo-scoped Codex skill under
`.agents/skills/bulk-atacseq-agent/` is an optional adapter. The files in this
directory remain the controlling provider-neutral contract.

## Supported modes

| Mode | Intended use | Computation allowed |
|---|---|---|
| `case_study_audit` | Inspect the published internship example | No rerun by default |
| `new_study_dry_run` | Validate a new study and prepare a staged plan | Configuration checks only |
| `new_study_execution` | Advance one approved stage | Only the explicitly approved stage |

## Recommended lifecycle

1. Copy `config/templates/samples.template.tsv` and
   `config/templates/analysis_parameters.agent.template.json` to new files
   under the Git-ignored `config/local/` directory.
2. Copy `analysis_status.template.json` to the Git-ignored
   `agent/local_status/` directory.
3. Replace every placeholder and add the actual study metadata.
4. Start the agent in `new_study_dry_run` mode using
   [`example_request.md`](example_request.md).
5. Run the configuration validator. Resolve all errors and review every warning.
6. Let the agent identify missing human decisions and propose a stage plan.
7. Approve at most one computational stage at a time.
8. After execution, require the stage's outputs and validation evidence before
   approving its successor.
9. Preserve the status record, parameters, commands, versions, summaries,
   limitations, and data-access instructions.

## Validation command

From the repository root:

```bash
python scripts/validate_agent_configuration.py \
  --samples config/local/my_samples.tsv \
  --parameters config/local/my_analysis_parameters.json \
  --status agent/local_status/my_analysis_status.json
```

Add `--check-files` only in an environment where the recorded FASTQ paths are
expected to exist. Add `--strict` when warnings should also produce a nonzero
exit status.

The validator does not prove that the experiment is scientifically adequate.
It checks structural invariants and common inconsistencies so human review can
focus on the design and evidence.

## Status vocabulary

- `not_started`: no work has been performed.
- `awaiting_human`: a required scientific choice or authorization is missing.
- `approved`: the stage is authorized but has not begun.
- `running`: the authorized stage is currently executing.
- `validation_failed`: execution finished but required checks did not pass.
- `complete`: outputs and validation evidence are recorded.
- `blocked`: progress requires an external input or infrastructure change.
- `not_applicable`: an optional stage was deliberately skipped with a reason.

## Portability

Agents that do not recognize `AGENTS.md` or Codex skills can still be instructed
to read this file and `workflow_contract.yaml`. JSON, JSON Schema, TSV, Markdown,
and the validator are intentionally independent of a particular model provider.

For Codex, open the repository as the working directory and invoke
`$bulk-atacseq-agent`, or describe a matching ATAC-seq task and allow normal
skill discovery. Codex documents repository-root
[`AGENTS.md`](https://learn.chatgpt.com/docs/agent-configuration/agents-md)
instructions and repo-scoped skills under `.agents/skills` in its
[skill guidance](https://learn.chatgpt.com/docs/build-skills). Other agents
should be told explicitly to read the root instructions and this directory.

This interface is designed for conventional bulk, paired-end ATAC-seq. A user
must choose a different validated workflow for single-cell or multiome data.
