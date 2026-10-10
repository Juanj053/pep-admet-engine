"""
Data models for peptide protease cleavage profiling.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional


@dataclass(frozen=True)
class CleavageSite:
    """Represents a detected proteolytic cleavage site."""
    site_index: int  # 1-based index (bond between site_index and site_index + 1)
    protease_name: str
    p1_residue: str
    p1_prime_residue: str
    cleavage_window: str  # e.g., P4P3P2P1|P1'P2'P3'P4'
    specificity_class: str = "canonical"
    cleavage_efficiency: float = 1.0

    def to_dict(self) -> Dict[str, object]:
        return {
            "site_index": self.site_index,
            "protease_name": self.protease_name,
            "p1_residue": self.p1_residue,
            "p1_prime_residue": self.p1_prime_residue,
            "cleavage_window": self.cleavage_window,
            "specificity_class": self.specificity_class,
            "cleavage_efficiency": self.cleavage_efficiency,
        }


@dataclass(frozen=True)
class PeptideFragment:
    """Represents a peptide fragment resulting from digestion."""
    fragment_index: int  # 1-based
    start_pos: int  # 1-based
    end_pos: int  # 1-based
    sequence: str
    length: int
    molecular_weight: float
    isoelectric_point: float

    def to_dict(self) -> Dict[str, object]:
        return {
            "fragment_index": self.fragment_index,
            "start_pos": self.start_pos,
            "end_pos": self.end_pos,
            "sequence": self.sequence,
            "length": self.length,
            "molecular_weight": round(self.molecular_weight, 2),
            "isoelectric_point": round(self.isoelectric_point, 2),
        }


@dataclass(frozen=True)
class ProteaseSpec:
    """Enzyme cleavage specificity and biological metadata."""
    name: str
    family: str  # Serine, Aspartic, Cysteine, Metalloprotease, Dipeptidyl
    description: str
    biological_context: str  # serum, digestive, lysosomal, plasma
    optimum_ph: float
    p1_allowed: List[str] = field(default_factory=list)
    p1_disallowed: List[str] = field(default_factory=list)
    p1_prime_disallowed: List[str] = field(default_factory=lambda: ["P"])
    p1_prime_allowed: List[str] = field(default_factory=list)
    n_terminal_only: bool = False
    c_terminal_only: bool = False


@dataclass
class CleavageResult:
    """Comprehensive cleavage report across selected proteases."""
    peptide_id: str
    sequence: str
    length: int
    molecular_weight: float
    proteases_analyzed: List[str]
    cleavage_sites: List[CleavageSite] = field(default_factory=list)
    fragments: List[PeptideFragment] = field(default_factory=list)
    cleavage_count: int = 0
    proteolytic_vulnerability_index: float = 0.0  # 0 to 100
    proteolytic_resistance_score: float = 100.0  # 0 to 100
    hotspot_residues: List[int] = field(default_factory=list)
    cleavage_map: Dict[int, List[str]] = field(default_factory=dict)

    def to_dict(self) -> Dict[str, object]:
        return {
            "peptide_id": self.peptide_id,
            "sequence": self.sequence,
            "length": self.length,
            "molecular_weight": round(self.molecular_weight, 2),
            "proteases_analyzed": self.proteases_analyzed,
            "cleavage_count": self.cleavage_count,
            "proteolytic_vulnerability_index": round(self.proteolytic_vulnerability_index, 2),
            "proteolytic_resistance_score": round(self.proteolytic_resistance_score, 2),
            "hotspot_residues": self.hotspot_residues,
            "cleavage_map": self.cleavage_map,
            "cleavage_sites": [s.to_dict() for s in self.cleavage_sites],
            "fragments": [f.to_dict() for f in self.fragments],
        }
