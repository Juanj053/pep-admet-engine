"""
PEP-ADMET Engine: Rational Stabilization Advisor
===============================================
Formulates medicinal chemistry and peptidomimetic engineering solutions
to overcome proteolytic vulnerability and pharmacokinetic liabilities.
"""

from __future__ import annotations

from typing import Dict, List, Optional

from pep_admet.models.cleavage import CleavageResult
from pep_admet.models.halflife import HalfLifePrediction
from pep_admet.models.permeability import PermeabilityMetrics
from pep_admet.models.stabilization import (
    StabilizationPlan,
    StabilizationRecommendation,
)


class StabilizationAdvisor:
    """Generates rational engineering recommendations for peptide drug leads."""

    def formulate_plan(
        self,
        cleavage: CleavageResult,
        halflife: Optional[HalfLifePrediction] = None,
        permeability: Optional[PermeabilityMetrics] = None,
    ) -> StabilizationPlan:
        seq = cleavage.sequence
        length = len(seq)
        recommendations: List[StabilizationRecommendation] = []
        engineered_variants: List[Dict[str, str]] = []

        # 1. N-terminal Capping Strategy
        # Check for DPP-4 or general aminopeptidase vulnerability
        dpp4_cut = any(s.protease_name.startswith("Dipeptidyl") for s in cleavage.cleavage_sites)
        if dpp4_cut:
            recommendations.append(
                StabilizationRecommendation(
                    strategy_type="N-terminal capping / substitution",
                    target_site="N-terminus (Residues 1-2)",
                    modification="N-terminal Acetylation (Ac-) or D-Ala / Aib at position 2",
                    rationale="Blocks DPP-IV dipeptidyl cleavage recognition while preserving N-terminal bioactivity.",
                    predicted_half_life_gain_fold=3.5,
                    potential_binding_impact="Minimal / Low risk",
                    priority=1,
                )
            )
            # generate variant
            mod_seq = f"Ac-{seq}"
            engineered_variants.append({"name": "N-term-Acetylated", "sequence": mod_seq})
        else:
            recommendations.append(
                StabilizationRecommendation(
                    strategy_type="N-terminal capping",
                    target_site="N-terminus",
                    modification="N-Acetylation (Ac-)",
                    rationale="Protects against circulating serum aminopeptidases and neutralizes positive charge.",
                    predicted_half_life_gain_fold=1.8,
                    potential_binding_impact="Minimal / Low risk",
                    priority=2,
                )
            )

        # 2. C-terminal Capping Strategy
        recommendations.append(
            StabilizationRecommendation(
                strategy_type="C-terminal capping",
                target_site="C-terminus",
                modification="C-terminal Amidation (-NH2)",
                rationale="Eliminates negative carboxylate, blocks carboxypeptidases, mimics native peptide hormones.",
                predicted_half_life_gain_fold=2.2,
                potential_binding_impact="Minimal / Low risk",
                priority=1,
            )
        )

        # 3. Scissile Bond Hardening (D-amino acids or N-methylation)
        if cleavage.hotspot_residues:
            for hp in cleavage.hotspot_residues[:3]:
                res = seq[hp - 1]
                recommendations.append(
                    StabilizationRecommendation(
                        strategy_type="Scissile bond stereoinversion / N-methylation",
                        target_site=f"Residue {hp} ({res})",
                        modification=f"D-{res} substitution or N-alpha-methyl-{res}",
                        rationale=f"Inverts scissile bond stereochemistry at high-frequency cleavage site {hp}, abolishing protease active-site fit.",
                        predicted_half_life_gain_fold=4.0,
                        potential_binding_impact="Moderate / Verify epitope",
                        priority=1,
                    )
                )

        # 4. Cyclization Engineering (Disulfide, Stapling, or Backbone)
        if length >= 6 and length <= 25:
            recommendations.append(
                StabilizationRecommendation(
                    strategy_type="Cyclization engineering",
                    target_site="Termini or i, i+4 flank",
                    modification="Cys-Cys disulfide bridging or hydrocarbon stapling",
                    rationale="Enforces bioactive conformation, shields backbone amides from endopeptidases, reduces polar surface area.",
                    predicted_half_life_gain_fold=5.5,
                    potential_binding_impact="Moderate / Verify epitope",
                    priority=2,
                )
            )

        # 5. Serum Albumin (HSA) Binding / Lipidation (Semaglutide design pattern)
        if halflife and halflife.estimated_t_half_minutes < 120:
            recommendations.append(
                StabilizationRecommendation(
                    strategy_type="PK Extension / Lipidation",
                    target_site="Lys side-chain or C-terminus",
                    modification="C16/C18 fatty diacid (palmitoyl) conjugation via mini-PEG/gamma-Glu spacer",
                    rationale="Promotes high-affinity reversible binding to Human Serum Albumin, surpassing glomerular renal filtration threshold.",
                    predicted_half_life_gain_fold=15.0,
                    potential_binding_impact="Low risk if conjugated outside core pharmacophore",
                    priority=1,
                )
            )

        summary = (
            f"Analyzed {length}-mer peptide '{cleavage.peptide_id}'. Identified {cleavage.cleavage_count} cleavage sites "
            f"with {len(cleavage.hotspot_residues)} proteolytic hotspots. "
            f"Top engineering priorities include terminal capping, scissile bond chiral inversion, and albumin-binding lipidation."
        )

        return StabilizationPlan(
            peptide_id=cleavage.peptide_id,
            wildtype_sequence=seq,
            vulnerability_summary=summary,
            recommendations=recommendations,
            engineered_variants=engineered_variants,
        )
