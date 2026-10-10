#!/usr/bin/env bash
set -euo pipefail

echo "=================================================="
echo "Running PEP-ADMET Engine Hermetic Test Harness (CI)"
echo "=================================================="

# Verify config schema validation
python3 -c "
import yaml, jsonschema
from pathlib import Path

config = yaml.safe_load(Path('.test/config/config.yaml').read_text())
schema = yaml.safe_load(Path('config/schemas/config.schema.yaml').read_text())
jsonschema.validate(instance=config, schema=schema)
print('Config schema validation: PASS')
"

# Run Snakemake dry-run
snakemake --configfile .test/config/config.yaml -n -q

echo "All CI test harness steps passed successfully."
