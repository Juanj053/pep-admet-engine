"""
Data models for membrane permeability and oral bioavailability.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional


@dataclass
class PermeabilityMetrics:
    """Physicochemical and transcellular permeability descriptors."""
    peptide_id: str
    sequence: str
    molecular_weight: float
    tpsa: float  # Topological Polar Surface Area (Å²)
    log_p: float  # Lipophilicity index (GRAVY / Wimley-White)
    hbd_count: int  # Hydrogen Bond Donors
    hba_count: int  # Hydrogen Bond Acceptors
    rotatable_bonds_count: int
    net_charge_ph74: float
    caco2_log_papp_cm_s: float  # Predicted Caco-2 permeability
    pampa_permeability_class: str  # "Low", "Medium", "High"
    blood_brain_barrier_log_bb: float  # Predicted BBB penetration
    cell_penetrating_score: float  # 0.0 to 1.0 (CPP propensity)
    oral_bioavailability_score: float  # 0.0 to 100.0 (bRo5 / Veber compliance)

    def to_dict(self) -> Dict[str, object]:
        return {
            "peptide_id": self.peptide_id,
            "sequence": self.sequence,
            "molecular_weight": round(self.molecular_weight, 2),
            "tpsa": round(self.tpsa, 2),
            "log_p": round(self.log_p, 2),
            "hbd_count": self.hbd_count,
            "hba_count": self.hba_count,
            "rotatable_bonds_count": self.rotatable_bonds_count,
            "net_charge_ph74": round(self.net_charge_ph74, 2),
            "caco2_log_papp_cm_s": round(self.caco2_log_papp_cm_s, 2),
            "pampa_permeability_class": self.pampa_permeability_class,
            "blood_brain_barrier_log_bb": round(self.blood_brain_barrier_log_bb, 2),
            "cell_penetrating_score": round(self.cell_penetrating_score, 2),
            "oral_bioavailability_score": round(self.oral_bioavailability_score, 2),
        }
