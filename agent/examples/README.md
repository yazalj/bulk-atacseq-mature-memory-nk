# Configuration-only example

These files demonstrate a structurally valid four-sample paired design. They do
not identify a real experiment and no sequencing data are included.

From the repository root, validate them with:

```bash
python scripts/validate_agent_configuration.py \
  --samples agent/examples/minimal_samples.tsv \
  --parameters agent/examples/minimal_parameters.json \
  --status agent/examples/minimal_status.json
```

The example remains in dry-run mode with human decisions pending. That is
expected; configuration validity is not scientific approval.
