"""
Data models for rational stabilization and peptidomimetic engineering.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass
class StabilizationRecommendation:
    """Specific chemical modification advice to improve in-vivo stability."""
    strategy_type: str  # "N-terminal capping", "C-terminal capping", "D-amino acid substitution", "Cyclization", "Lipidation", "Backbone modification"
    target_site: str  # "N-terminus", "C-terminus", "Residue 5 (Lys)", "Disulfide loop"
    modification: str  # "N-Acetylation (Ac-)", "C-Amidation (-NH2)", "D-Phe at P1", "C16 Palmitoyl"
    rationale: str  # Why this solves degradation or PK failure
    predicted_half_life_gain_fold: float
    potential_binding_impact: str  # "Minimal / Low risk", "Moderate / Verify epitope", "Substantial"
    priority: int  # 1 (Highest), 2, 3

    def to_dict(self) -> Dict[str, object]:
        return {
            "strategy_type": self.strategy_type,
            "target_site": self.target_site,
            "modification": self.modification,
            "rationale": self.rationale,
            "predicted_half_life_gain_fold": round(self.predicted_half_life_gain_fold, 2),
            "potential_binding_impact": self.potential_binding_impact,
            "priority": self.priority,
        }


@dataclass
class StabilizationPlan:
    """Full stabilization roadmap for a candidate peptide."""
    peptide_id: str
    wildtype_sequence: str
    vulnerability_summary: str
    recommendations: List[StabilizationRecommendation]
    engineered_variants: List[Dict[str, str]]  # list of modified analogs

    def to_dict(self) -> Dict[str, object]:
        return {
            "peptide_id": self.peptide_id,
            "wildtype_sequence": self.wildtype_sequence,
            "vulnerability_summary": self.vulnerability_summary,
            "recommendations": [r.to_dict() for r in self.recommendations],
            "engineered_variants": self.engineered_variants,
        }
