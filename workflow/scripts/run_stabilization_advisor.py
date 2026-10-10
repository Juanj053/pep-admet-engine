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
logger = logging.getLogger("stabilization_advisor")

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

from pep_admet.core.stabilization_advisor import StabilizationAdvisor
from pep_admet.models.cleavage import CleavageResult, CleavageSite, PeptideFragment
from pep_admet.models.halflife import HalfLifePrediction
from pep_admet.models.permeability import PermeabilityMetrics

sample = snakemake.wildcards.sample

with open(snakemake.input.cleavage_json, "r", encoding="utf-8") as f:
    c_raw = json.load(f)
with open(snakemake.input.halflife_json, "r", encoding="utf-8") as f:
    h_raw = json.load(f)
with open(snakemake.input.permeability_json, "r", encoding="utf-8") as f:
    p_raw = json.load(f)

cleavage = CleavageResult(
    peptide_id=c_raw["peptide_id"],
    sequence=c_raw["sequence"],
    length=c_raw["length"],
    molecular_weight=c_raw["molecular_weight"],
    proteases_analyzed=c_raw["proteases_analyzed"],
    cleavage_sites=[CleavageSite(**s) for s in c_raw["cleavage_sites"]],
    fragments=[PeptideFragment(**fg) for fg in c_raw["fragments"]],
    cleavage_count=c_raw["cleavage_count"],
    proteolytic_vulnerability_index=c_raw["proteolytic_vulnerability_index"],
    proteolytic_resistance_score=c_raw["proteolytic_resistance_score"],
    hotspot_residues=c_raw["hotspot_residues"],
    cleavage_map=c_raw["cleavage_map"],
)

halflife = HalfLifePrediction(**h_raw)
permeability = PermeabilityMetrics(**p_raw)

logger.info(f"Formulating stabilization plan for {sample}...")
advisor = StabilizationAdvisor()
plan = advisor.formulate_plan(cleavage, halflife, permeability)

out_json = Path(snakemake.output.stabilization_json)
out_json.parent.mkdir(parents=True, exist_ok=True)
with open(out_json, "w", encoding="utf-8") as f:
    json.dump(plan.to_dict(), f, indent=2)

logger.info(f"Stabilization roadmap formulated with {len(plan.recommendations)} actionable strategies.")
