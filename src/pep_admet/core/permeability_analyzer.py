"""
PEP-ADMET Engine: Membrane Permeability & Bioavailability Analyzer
==================================================================
Computes topological polar surface area (TPSA), logP, HBD/HBA,
predicted Caco-2 log Papp, Blood-Brain Barrier (log BB), and CPP score.
"""

from __future__ import annotations

import math
from typing import Dict

from pep_admet.models.permeability import PermeabilityMetrics

# TPSA contributions per residue (Å²)
RESIDUE_TPSA: Dict[str, float] = {
    "A": 29.1, "R": 73.2, "N": 63.3, "D": 66.4, "C": 29.1,
    "E": 66.4, "Q": 63.3, "G": 29.1, "H": 42.0, "I": 29.1,
    "L": 29.1, "K": 55.1, "M": 29.1, "F": 29.1, "P": 20.3,
    "S": 49.3, "T": 49.3, "W": 44.9, "Y": 49.3, "V": 29.1,
}

# Kyte-Doolittle Hydropathy values
HYDROPATHY: Dict[str, float] = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
    "E": -3.5, "Q": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
    "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
    "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}

# Hydrogen bond donor / acceptor counts per residue
HBD_MAP: Dict[str, int] = {
    "A": 1, "R": 4, "N": 2, "D": 1, "C": 1, "E": 1, "Q": 2,
    "G": 1, "H": 2, "I": 1, "L": 1, "K": 2, "M": 1, "F": 1,
    "P": 0, "S": 2, "T": 2, "W": 2, "Y": 2, "V": 1,
}
HBA_MAP: Dict[str, int] = {
    "A": 1, "R": 2, "N": 2, "D": 2, "C": 1, "E": 2, "Q": 2,
    "G": 1, "H": 2, "I": 1, "L": 1, "K": 1, "M": 1, "F": 1,
    "P": 1, "S": 2, "T": 2, "W": 1, "Y": 2, "V": 1,
}


class MembranePermeabilityAnalyzer:
    """Analyzes passive transcellular permeability and physicochemical drug-likeness."""

    def analyze(self, sequence: str, peptide_id: str = "peptide_candidate") -> PermeabilityMetrics:
        clean_seq = sequence.upper().strip()
        length = max(1, len(clean_seq))

        # Molecular weight
        mw = sum(RESIDUE_TPSA.get(aa, 30.0) * 2.5 for aa in clean_seq)  # approx scale
        from pep_admet.core.protease_engine import AA_RESIDUE_MASSES, WATER_MASS
        mw = sum(AA_RESIDUE_MASSES.get(aa, 110.0) for aa in clean_seq) + WATER_MASS

        # TPSA calculation (sum + terminal groups ~37.3 Å²)
        tpsa = sum(RESIDUE_TPSA.get(aa, 30.0) for aa in clean_seq) + 37.3

        # GRAVY / logP index
        gravy = sum(HYDROPATHY.get(aa, 0.0) for aa in clean_seq) / length

        # HBD & HBA counts
        hbd = sum(HBD_MAP.get(aa, 1) for aa in clean_seq) + 2  # + terminal N-term & C-term
        hba = sum(HBA_MAP.get(aa, 1) for aa in clean_seq) + 2

        # Rotatable bonds (~ 2 per standard residue + sidechain rotatable bonds)
        rot_bonds = length * 2 + sum(1 for aa in clean_seq if aa in ["K", "R", "Q", "E", "M"])

        # Net charge at pH 7.4
        pos_charge = clean_seq.count("K") + clean_seq.count("R") + (0.1 * clean_seq.count("H")) + 1.0  # N-term
        neg_charge = clean_seq.count("D") + clean_seq.count("E") + 1.0  # C-term
        net_charge = pos_charge - neg_charge

        # Cell Penetrating Peptide (CPP) score: high Arg/Lys density
        cationic_density = (clean_seq.count("R") * 1.5 + clean_seq.count("K")) / length
        cpp_score = min(1.0, cationic_density * 2.0)

        # Caco-2 Log Papp estimation (empirical QSAR: penalized by TPSA, HBD, MW)
        # Log Papp values typically range between -7.5 (low) to -5.0 (high) cm/s
        caco2_log_papp = -5.0 - (0.005 * tpsa) - (0.12 * hbd) + (0.35 * gravy)
        caco2_log_papp = max(-8.5, min(-4.5, caco2_log_papp))

        # PAMPA class
        if caco2_log_papp > -5.8:
            pampa_class = "High"
        elif caco2_log_papp > -6.8:
            pampa_class = "Medium"
        else:
            pampa_class = "Low"

        # Blood-Brain Barrier (log BB)
        # log BB = 0.152 logP - 0.0148 TPSA + 0.139
        log_bb = (0.152 * gravy) - (0.0148 * tpsa) + 0.139
        log_bb = max(-3.0, min(1.0, log_bb))

        # Oral Bioavailability score (0 to 100): bRo5 rules
        # Penalties for high TPSA, excessive HBD, excessive rotatable bonds
        oral_score = 100.0 - (tpsa * 0.15) - (hbd * 2.5) - (rot_bonds * 0.8) + (gravy * 5.0)
        oral_score = max(0.0, min(100.0, oral_score))

        return PermeabilityMetrics(
            peptide_id=peptide_id,
            sequence=clean_seq,
            molecular_weight=mw,
            tpsa=tpsa,
            log_p=gravy,
            hbd_count=hbd,
            hba_count=hba,
            rotatable_bonds_count=rot_bonds,
            net_charge_ph74=net_charge,
            caco2_log_papp_cm_s=caco2_log_papp,
            pampa_permeability_class=pampa_class,
            blood_brain_barrier_log_bb=log_bb,
            cell_penetrating_score=cpp_score,
            oral_bioavailability_score=oral_score,
        )
