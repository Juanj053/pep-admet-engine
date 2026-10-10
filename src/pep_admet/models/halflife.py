"""
Data models for peptide plasma half-life and pharmacokinetics.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass
class HalfLifePrediction:
    """Predicted systemic circulation half-life and clearance parameters."""
    peptide_id: str
    sequence: str
    length: int
    molecular_weight: float
    estimated_t_half_minutes: float
    estimated_t_half_category: str  # "ultra-short (<15 min)", "short (15-60 min)", "moderate (1-6 h)", "extended (>6 h)"
    renal_filtration_risk: str  # "High (<5 kDa)", "Moderate (5-30 kDa)", "Low (>30 kDa)"
    proteolytic_instability_penalty: float
    albumin_binding_propensity: float  # 0.0 to 1.0
    predicted_clearance_ml_min_kg: float
    predicted_vss_l_kg: float

    def to_dict(self) -> Dict[str, object]:
        return {
            "peptide_id": self.peptide_id,
            "sequence": self.sequence,
            "length": self.length,
            "molecular_weight": round(self.molecular_weight, 2),
            "estimated_t_half_minutes": round(self.estimated_t_half_minutes, 1),
            "estimated_t_half_category": self.estimated_t_half_category,
            "renal_filtration_risk": self.renal_filtration_risk,
            "proteolytic_instability_penalty": round(self.proteolytic_instability_penalty, 2),
            "albumin_binding_propensity": round(self.albumin_binding_propensity, 2),
            "predicted_clearance_ml_min_kg": round(self.predicted_clearance_ml_min_kg, 2),
            "predicted_vss_l_kg": round(self.predicted_vss_l_kg, 2),
        }
