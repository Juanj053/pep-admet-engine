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
logger = logging.getLogger("reporting")

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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pep_admet.core.report_builder import ADMETReportBuilder
from pep_admet.models.cleavage import CleavageResult, CleavageSite, PeptideFragment
from pep_admet.models.halflife import HalfLifePrediction
from pep_admet.models.permeability import PermeabilityMetrics
from pep_admet.models.stabilization import StabilizationPlan, StabilizationRecommendation

samples = list(snakemake.config["peptides"].keys())
builder = ADMETReportBuilder()
profiles = []

for sample in samples:
    with open(f"results/{sample}/cleavage_profile.json", "r", encoding="utf-8") as f:
        c_raw = json.load(f)
    with open(f"results/{sample}/plasma_halflife.json", "r", encoding="utf-8") as f:
        h_raw = json.load(f)
    with open(f"results/{sample}/membrane_permeability.json", "r", encoding="utf-8") as f:
        p_raw = json.load(f)
    with open(f"results/{sample}/stabilization_plan.json", "r", encoding="utf-8") as f:
        s_raw = json.load(f)

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
    stabilization = StabilizationPlan(
        peptide_id=s_raw["peptide_id"],
        wildtype_sequence=s_raw["wildtype_sequence"],
        vulnerability_summary=s_raw["vulnerability_summary"],
        recommendations=[StabilizationRecommendation(**r) for r in s_raw["recommendations"]],
        engineered_variants=s_raw["engineered_variants"],
    )

    prof = builder.build_profile(cleavage, halflife, permeability, stabilization)
    profiles.append(prof)

    # Export per-sample JSON
    builder.export_json(prof, f"results/{sample}/admet_profile.json")

    # Generate publication-grade 4-panel visual analytics figure
    if snakemake.params.generate_plots:
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        fig.suptitle(f"PEP-ADMET Analytics: {sample} (Score: {prof.overall_admet_score:.1f}/100)", fontsize=16, fontweight="bold")

        # Panel A: Cleavage site distribution across peptide sequence
        ax_a = axes[0, 0]
        positions = list(range(1, prof.cleavage.length))
        cuts_per_pos = [len(prof.cleavage.cleavage_map.get(p, [])) for p in positions]
        ax_a.bar(positions, cuts_per_pos, color="#e74c3c", alpha=0.85, edgecolor="black")
        ax_a.set_title("Proteolytic Cleavage Frequency per Scissile Bond", fontweight="bold")
        ax_a.set_xlabel("Residue Position (P1)")
        ax_a.set_ylabel("Number of Cleaving Proteases")
        ax_a.grid(True, linestyle="--", alpha=0.4)

        # Panel B: Simulated Degradation Kinetics C(t) = C0 * exp(-k*t)
        ax_b = axes[0, 1]
        t = np.linspace(0, 120, 200)
        t_half = max(5.0, prof.halflife.estimated_t_half_minutes)
        k_elim = np.log(2) / t_half
        c_wt = 100.0 * np.exp(-k_elim * t)
        # Stabilized variant simulation (3x - 5x t1/2 gain)
        k_stab = k_elim / 4.0
        c_stab = 100.0 * np.exp(-k_stab * t)
        ax_b.plot(t, c_wt, "r-", label=f"Wild Type (t1/2 ~ {t_half:.0f}m)", linewidth=2)
        ax_b.plot(t, c_stab, "g--", label=f"Stabilized Analog (t1/2 ~ {t_half*4:.0f}m)", linewidth=2)
        ax_b.set_title("In-Silico Systemic Elimination Kinetics", fontweight="bold")
        ax_b.set_xlabel("Time (minutes)")
        ax_b.set_ylabel("Intact Peptide Concentration (%)")
        ax_b.legend(loc="upper right")
        ax_b.grid(True, linestyle="--", alpha=0.4)

        # Panel C: ADMET Radar Chart (normalized properties)
        ax_c = axes[1, 0]
        categories = ["Stability", "Half-Life", "Permeability", "Oral", "Albumin Aff"]
        values = [
            prof.cleavage.proteolytic_resistance_score,
            min(100.0, (prof.halflife.estimated_t_half_minutes / 180.0) * 100.0),
            max(0.0, min(100.0, (prof.permeability.caco2_log_papp_cm_s + 8.0) * 25.0)),
            prof.permeability.oral_bioavailability_score,
            prof.halflife.albumin_binding_propensity * 100.0,
        ]
        bars = ax_c.barh(categories, values, color="#3498db", alpha=0.8, edgecolor="black")
        ax_c.set_xlim(0, 100)
        ax_c.set_title("ADMET Profile Pillars (0 - 100)", fontweight="bold")
        ax_c.grid(True, linestyle="--", alpha=0.4)
        for bar in bars:
            ax_c.text(bar.get_width() + 1.5, bar.get_y() + bar.get_height()/2, f"{bar.get_width():.1f}", va="center", fontweight="bold")

        # Panel D: Stabilization Strategies Summary
        ax_d = axes[1, 1]
        ax_d.axis("off")
        recs_lines = ["Top Stabilization Recommendations:", ""]
        for idx, r in enumerate(prof.stabilization.recommendations[:4], start=1):
            recs_lines.append(f"{idx}. [{r.strategy_type}]")
            recs_lines.append(f"   Mod: {r.modification}")
            recs_lines.append(f"   Gain: {r.predicted_half_life_gain_fold:.1f}x t1/2 | Risk: {r.potential_binding_impact}")
            recs_lines.append("")
        recs_text = chr(10).join(recs_lines)
        ax_d.text(0.05, 0.95, recs_text, transform=ax_d.transAxes, verticalalignment="top", fontfamily="monospace", fontsize=9, bbox=dict(boxstyle="round,pad=0.5", facecolor="#f8f9fa", edgecolor="#bdc3c7"))

        plt.tight_layout()
        plt.savefig(f"results/{sample}/admet_analytics.png", dpi=snakemake.params.dpi)
        plt.close()

# Export multi-sample summary CSV
builder.export_summary_csv(profiles, snakemake.output.summary_csv)

# Export Global JSON
global_data = {
    "samples_count": len(profiles),
    "profiles": [p.to_dict() for p in profiles],
}
with open(snakemake.output.global_report_json, "w", encoding="utf-8") as f:
    json.dump(global_data, f, indent=2)

logger.info(f"Reporting completed successfully across {len(profiles)} peptides.")
