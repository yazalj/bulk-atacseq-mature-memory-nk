# Configuration templates for a new study

Copy these files into a new project-specific configuration directory; do not
edit the templates in place and assume they are analysis-ready.

- `samples.template.tsv` records biological identity, experimental factors,
  sequencing layout, and paired read paths.
- `input_manifest.template.tsv` records input roles, provenance, and checksums.
- `analysis_parameters.template.yaml` makes study-specific choices explicit.
- `analysis_parameters.agent.template.json` provides the same role in a format
  that can be checked deterministically by the AI-agent validator.

Replace every example row, placeholder path, `REVIEW_ME`, and
`REVIEW_FROM_*` value. Add columns needed for the actual design, such as sex,
time point, tissue, stimulation, or sequencing batch. Remove a column only when
it is genuinely irrelevant, not because its value is unknown.

For an AI-assisted local run, place filled copies under `config/local/`. That
directory is ignored by Git so private absolute paths are not published by
accident. The templates and configuration-only examples remain tracked.

These general templates do not have the locked schema required by the
internship case-study scripts under `workflow/`. If adapting those scripts,
revise their manifest parsing, expected dimensions, factor levels, and
validation assertions together with the new study design.
