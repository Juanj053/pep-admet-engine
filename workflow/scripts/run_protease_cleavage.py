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
logger = logging.getLogger("protease_cleavage")

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

from pep_admet.core.protease_engine import ProteaseEngine

sample = snakemake.wildcards.sample
seq = snakemake.params.sequence
proteases = snakemake.params.proteases

logger.info(f"Running protease cleavage analysis for {sample} ({len(seq)} AA)...")
engine = ProteaseEngine()
result = engine.analyze(sequence=seq, peptide_id=sample, proteases=proteases)

# Save JSON
out_json = Path(snakemake.output.cleavage_json)
out_json.parent.mkdir(parents=True, exist_ok=True)
with open(out_json, "w", encoding="utf-8") as f:
    json.dump(result.to_dict(), f, indent=2)

# Save Fragments CSV
out_csv = Path(snakemake.output.fragments_csv)
fieldnames = ["fragment_index", "start_pos", "end_pos", "sequence", "length", "molecular_weight", "isoelectric_point"]
with open(out_csv, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for frag in result.fragments:
        writer.writerow(frag.to_dict())

logger.info(f"Protease analysis complete. Found {result.cleavage_count} cuts across {len(result.fragments)} fragments.")
