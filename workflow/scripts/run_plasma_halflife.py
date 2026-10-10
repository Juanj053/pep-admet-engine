import csv
import json
import logging
import sys
from pathlib import Path

if hasattr(snakemake, "log") and snakemake.log:
    log_path = Path(snakemake.log[0])
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_file = open(log_path, "w", encoding="utf-8")
    sys.stdout = log_file
    sys.stderr = log_file

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("plasma_halflife")

def _resolve_repo_root() -> Path:
    snk = globals().get("snakemake", None)
    if snk is not None and hasattr(snk, "params"):
        wb = getattr(snk.params, "workflow_basedir", None)
        if wb:
            cand = Path(wb).resolve()
            for p in [cand, cand.parent] + list(cand.parents):
                if (p / "src" / "pep_admet").is_dir():
                    return p
    script_file = globals().get("__file__", None)
    if script_file:
        cand = Path(script_file).resolve()
        for p in [cand.parent] + list(cand.parents):
            if (p / "src" / "pep_admet").is_dir():
                return p
    cand = Path.cwd().resolve()
    for p in [cand] + list(cand.parents):
        if (p / "src" / "pep_admet").is_dir():
            return p
    raise RuntimeError("Unable to locate repository root containing src/pep_admet")

repo_root = _resolve_repo_root()
if str(repo_root / "src") not in sys.path:
    sys.path.insert(0, str(repo_root / "src"))

from pep_admet.core.halflife_predictor import PlasmaHalfLifePredictor

sample = snakemake.wildcards.sample
seq = snakemake.params.sequence

# Read cleavage result to extract proteolytic vulnerability index
in_cleavage = Path(snakemake.input.cleavage_json)
with open(in_cleavage, "r", encoding="utf-8") as f:
    cleavage_data = json.load(f)
pvi = cleavage_data.get("proteolytic_vulnerability_index", 20.0)

logger.info(f"Predicting plasma half-life for {sample} (PVI={pvi:.1f})...")
predictor = PlasmaHalfLifePredictor()
prediction = predictor.predict(sequence=seq, peptide_id=sample, proteolytic_vulnerability_index=pvi)

out_json = Path(snakemake.output.halflife_json)
out_json.parent.mkdir(parents=True, exist_ok=True)
with open(out_json, "w", encoding="utf-8") as f:
    json.dump(prediction.to_dict(), f, indent=2)

logger.info(f"Predicted half-life: {prediction.estimated_t_half_minutes:.1f} min ({prediction.estimated_t_half_category}).")
