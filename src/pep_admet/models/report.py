"""
Composite data models for end-to-end ADMET reports.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

from .cleavage import CleavageResult
from .halflife import HalfLifePrediction
from .permeability import PermeabilityMetrics
from .stabilization import StabilizationPlan


@dataclass
class PeptideADMETProfile:
    """Integrated ADMET profile combining all sub-engines."""
    peptide_id: str
    sequence: str
    cleavage: CleavageResult
    halflife: HalfLifePrediction
    permeability: PermeabilityMetrics
    stabilization: StabilizationPlan
    overall_admet_score: float  # 0 to 100

    def to_dict(self) -> Dict[str, object]:
        return {
            "peptide_id": self.peptide_id,
            "sequence": self.sequence,
            "overall_admet_score": round(self.overall_admet_score, 1),
            "cleavage": self.cleavage.to_dict(),
            "halflife": self.halflife.to_dict(),
            "permeability": self.permeability.to_dict(),
            "stabilization": self.stabilization.to_dict(),
        }
