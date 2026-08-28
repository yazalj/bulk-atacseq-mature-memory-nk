# Example request for an AI agent

Copy and adapt the text below when starting a new analysis.

> Read `AGENTS.md`, `agent/README.md`, `agent/workflow_contract.yaml`, and the
> relevant tutorial chapters. Work in `new_study_dry_run` mode. My sample
> registry is `config/local/my_samples.tsv`, my parameters are
> `config/local/my_analysis_parameters.json`, and my status file is
> `agent/local_status/my_analysis_status.json`. Validate the configuration without modifying
> raw or reference inputs and without running biological computation. Report
> structural errors, warnings, unresolved scientific decisions, and the next
> single stage that could be approved. Do not guess missing values or advance
> past a human-approval gate.

After reviewing the dry run, authorize one stage explicitly and name its output
directory. An instruction such as “continue the whole workflow” should not be
treated as approval to bypass unresolved decisions or failed validation gates.
