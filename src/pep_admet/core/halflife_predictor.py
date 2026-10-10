"""
PEP-ADMET Engine: Plasma Half-Life & Pharmacokinetics Predictor
==============================================================
Estimates in-vivo systemic clearance, volume of distribution,
renal filtration cutoff, and elimination half-life (t1/2).
"""

from __future__ import annotations

import math
from typing import Dict, Optional

from pep_admet.models.halflife import HalfLifePrediction

# Standard molecular weights
AA_MASSES: Dict[str, float] = {
    "A": 71.08, "R": 156.19, "N": 114.10, "D": 115.09, "C": 103.14,
    "E": 129.12, "Q": 128.13, "G": 57.05, "H": 137.14, "I": 113.16,
    "L": 113.16, "K": 128.17, "M": 131.19, "F": 147.18, "P": 97.12,
    "S": 87.08, "T": 101.11, "W": 186.21, "Y": 163.18, "V": 99.13,
}
WATER_MASS = 18.015


class PlasmaHalfLifePredictor:
    """Predicts plasma stability and systemic pharmacokinetic parameters."""

    def predict(
        self,
        sequence: str,
        peptide_id: str = "peptide_candidate",
        proteolytic_vulnerability_index: float = 20.0,
    ) -> HalfLifePrediction:
        clean_seq = sequence.upper().strip()
        length = len(clean_seq)
        mw = sum(AA_MASSES.get(aa, 110.0) for aa in clean_seq) + WATER_MASS

        # 1. Renal filtration threshold:
        # Glomerulus molecular weight pore cutoff ~ 30-50 kDa.
        # Below 5 kDa: rapid renal clearance (t1/2 ~ 5-20 min).
        # 5-30 kDa: intermediate clearance.
        # >30 kDa: minimal renal filtration.
        if mw < 5000:
            renal_risk = "High (<5 kDa)"
            base_renal_cl = 8.5  # mL/min/kg
        elif mw < 30000:
            renal_risk = "Moderate (5-30 kDa)"
            base_renal_cl = 3.2
        else:
            renal_risk = "Low (>30 kDa)"
            base_renal_cl = 0.5

        # 2. Proteolytic instability penalty
        # Higher PVI dramatically accelerates clearance via tissue/endothelial proteases
        proteolytic_penalty = proteolytic_vulnerability_index / 100.0  # 0 to 1
        proteolytic_cl = proteolytic_penalty * 15.0  # up to +15 mL/min/kg

        # 3. Serum Albumin Binding propensity
        # Hydrophobic patches (Trp, Phe, Leu clusters) increase non-covalent HSA affinity
        hydrophobic_count = sum(clean_seq.count(aa) for aa in ["W", "F", "L", "I"])
        fraction_hydrophobic = hydrophobic_count / max(1, length)
        # Albumin binding reduces free fraction, diminishing renal clearance
        albumin_binding = min(0.95, max(0.05, fraction_hydrophobic * 1.5))

        # Effective total clearance
        effective_cl = (base_renal_cl * (1.0 - albumin_binding * 0.7)) + proteolytic_cl
        effective_cl = max(0.1, effective_cl)

        # Steady-state volume of distribution (Vss in L/kg)
        # Small hydrophilic peptides distribute primarily into extracellular fluid (~0.2 - 0.4 L/kg)
        vss = 0.25 + (fraction_hydrophobic * 0.3)

        # Elimination half-life: t1/2 = ln(2) * Vss / CL (converted to minutes)
        # CL in mL/min/kg -> L/min/kg is CL / 1000
        cl_l_min = effective_cl / 1000.0
        t_half_minutes = (math.log(2) * vss) / cl_l_min

        # Half-life category
        if t_half_minutes < 15:
            cat = "ultra-short (<15 min)"
        elif t_half_minutes < 60:
            cat = "short (15-60 min)"
        elif t_half_minutes < 360:
            cat = "moderate (1-6 h)"
        else:
            cat = "extended (>6 h)"

        return HalfLifePrediction(
            peptide_id=peptide_id,
            sequence=clean_seq,
            length=length,
            molecular_weight=mw,
            estimated_t_half_minutes=t_half_minutes,
            estimated_t_half_category=cat,
            renal_filtration_risk=renal_risk,
            proteolytic_instability_penalty=proteolytic_penalty,
            albumin_binding_propensity=albumin_binding,
            predicted_clearance_ml_min_kg=effective_cl,
            predicted_vss_l_kg=vss,
        )
