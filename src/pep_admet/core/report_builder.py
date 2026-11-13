"""
PEP-ADMET Engine: Comprehensive ADMET Report Builder
===================================================
Aggregates individual analytics into integrated profiles, summary CSVs,
and formatted data structures for Snakemake HTML Reporting.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import TYPE_CHECKING

from pep_admet.models.cleavage import CleavageResult
from pep_admet.models.halflife import HalfLifePrediction
from pep_admet.models.permeability import PermeabilityMetrics
from pep_admet.models.report import PeptideADMETProfile
from pep_admet.models.stabilization import StabilizationPlan

if TYPE_CHECKING:
    from pep_admet.core.halflife_predictor import PlasmaHalfLifePredictor
    from pep_admet.core.permeability_analyzer import MembranePermeabilityAnalyzer
    from pep_admet.core.protease_engine import ProteaseEngine


class ADMETReportBuilder:
    """Integrator for multi-engine peptide ADMET analytics."""

    def build_profile(
        self,
        cleavage: CleavageResult,
        halflife: HalfLifePrediction,
        permeability: PermeabilityMetrics,
        stabilization: StabilizationPlan,
    ) -> PeptideADMETProfile:
        # Calculate composite ADMET score (0 to 100)
        # Weights: Stability (35%), Half-life (25%), Permeability/Oral (25%), Drug-likeness (15%)
        stability_comp = cleavage.proteolytic_resistance_score * 0.35
        halflife_comp = min(100.0, (halflife.estimated_t_half_minutes / 360.0) * 100.0) * 0.25
        perm_comp = permeability.oral_bioavailability_score * 0.25
        drug_comp = (100.0 - min(100.0, permeability.rotatable_bonds_count * 2.0)) * 0.15

        overall_score = max(0.0, min(100.0, stability_comp + halflife_comp + perm_comp + drug_comp))

        return PeptideADMETProfile(
            peptide_id=cleavage.peptide_id,
            sequence=cleavage.sequence,
            cleavage=cleavage,
            halflife=halflife,
            permeability=permeability,
            stabilization=stabilization,
            overall_admet_score=overall_score,
        )

    def export_json(self, profile: PeptideADMETProfile, output_path: str | Path) -> None:
        p = Path(output_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        with open(p, "w", encoding="utf-8") as f:
            json.dump(profile.to_dict(), f, indent=2)

    def compare_analogs(
        self,
        wildtype_profile: PeptideADMETProfile,
        protease_engine: ProteaseEngine,
        halflife_predictor: PlasmaHalfLifePredictor,
        permeability_analyzer: MembranePermeabilityAnalyzer,
    ) -> list[dict[str, float | str]]:
        """Screen and rank engineered peptidomimetic analogs against the wildtype profile."""
        comparisons: list[dict[str, float | str]] = []
        wt_t_half = max(1e-6, wildtype_profile.halflife.estimated_t_half_minutes)

        for entry in wildtype_profile.stabilization.engineered_variants:
            variant_id = entry.get("name", "variant")
            raw_seq = entry["sequence"]

            has_n_cap = raw_seq.startswith("Ac-")
            has_c_cap = raw_seq.endswith("-NH2")
            clean_seq = raw_seq
            if has_n_cap:
                clean_seq = clean_seq[len("Ac-") :]
            if has_c_cap:
                clean_seq = clean_seq[: -len("-NH2")]

            active_proteases = [
                p for p in protease_engine.available_proteases if not (has_n_cap and p == "dpp4")
            ]
            var_cleavage = protease_engine.analyze(
                clean_seq, peptide_id=variant_id, proteases=active_proteases
            )
            effective_pvi = max(
                2.0,
                var_cleavage.proteolytic_vulnerability_index
                * (0.75 if has_n_cap else 1.0)
                * (0.85 if has_c_cap else 1.0),
            )
            var_cleavage.proteolytic_vulnerability_index = effective_pvi
            var_cleavage.proteolytic_resistance_score = max(0.0, 100.0 - effective_pvi)

            var_halflife = halflife_predictor.predict(
                clean_seq,
                peptide_id=variant_id,
                proteolytic_vulnerability_index=effective_pvi,
            )
            var_perm = permeability_analyzer.analyze(clean_seq, peptide_id=variant_id)
            var_profile = self.build_profile(
                var_cleavage, var_halflife, var_perm, wildtype_profile.stabilization
            )
            comparisons.append(
                {
                    "variant_id": variant_id,
                    "sequence": raw_seq,
                    "proteolytic_vulnerability_index": round(effective_pvi, 2),
                    "estimated_t_half_minutes": round(var_halflife.estimated_t_half_minutes, 2),
                    "half_life_fold_vs_wt": round(
                        var_halflife.estimated_t_half_minutes / wt_t_half, 2
                    ),
                    "caco2_log_papp_cm_s": round(var_perm.caco2_log_papp_cm_s, 2),
                    "overall_admet_score": round(var_profile.overall_admet_score, 1),
                }
            )

        return comparisons

    def export_summary_csv(self, profiles: list[PeptideADMETProfile], output_path: str | Path) -> None:
        p = Path(output_path)
        p.parent.mkdir(parents=True, exist_ok=True)
        fieldnames = [
            "peptide_id",
            "sequence",
            "length",
            "molecular_weight",
            "overall_admet_score",
            "cleavage_count",
            "proteolytic_resistance_score",
            "estimated_t_half_minutes",
            "estimated_t_half_category",
            "renal_filtration_risk",
            "tpsa",
            "caco2_log_papp",
            "pampa_class",
            "oral_bioavailability_score",
            "top_stabilization_strategy",
        ]
        with open(p, "w", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for prof in profiles:
                top_strat = (
                    prof.stabilization.recommendations[0].modification
                    if prof.stabilization.recommendations
                    else "None"
                )
                writer.writerow({
                    "peptide_id": prof.peptide_id,
                    "sequence": prof.sequence,
                    "length": prof.cleavage.length,
                    "molecular_weight": round(prof.cleavage.molecular_weight, 2),
                    "overall_admet_score": round(prof.overall_admet_score, 1),
                    "cleavage_count": prof.cleavage.cleavage_count,
                    "proteolytic_resistance_score": round(prof.cleavage.proteolytic_resistance_score, 1),
                    "estimated_t_half_minutes": round(prof.halflife.estimated_t_half_minutes, 1),
                    "estimated_t_half_category": prof.halflife.estimated_t_half_category,
                    "renal_filtration_risk": prof.halflife.renal_filtration_risk,
                    "tpsa": round(prof.permeability.tpsa, 1),
                    "caco2_log_papp": round(prof.permeability.caco2_log_papp_cm_s, 2),
                    "pampa_class": prof.permeability.pampa_permeability_class,
                    "oral_bioavailability_score": round(prof.permeability.oral_bioavailability_score, 1),
                    "top_stabilization_strategy": top_strat,
                })
