"""
PEP-ADMET Engine: Protease Cleavage Profiler
============================================
Comprehensive in-silico enzymatic digestion engine for therapeutic peptides.
Implements recognition rules conforming to MEROPS classification, Schechter & Berger
subsite nomenclature, and empirical steric restriction filters.
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Set, Tuple

from pep_admet.models.cleavage import (
    CleavageResult,
    CleavageSite,
    PeptideFragment,
    ProteaseSpec,
)

# Standard average isotopic molecular masses for amino acid residues (minus H2O)
AA_RESIDUE_MASSES: Dict[str, float] = {
    "A": 71.0788,
    "R": 156.1875,
    "N": 114.1038,
    "D": 115.0886,
    "C": 103.1388,
    "E": 129.1155,
    "Q": 128.1307,
    "G": 57.0519,
    "H": 137.1411,
    "I": 113.1594,
    "L": 113.1594,
    "K": 128.1741,
    "M": 131.1926,
    "F": 147.1766,
    "P": 97.1167,
    "S": 87.0782,
    "T": 101.1051,
    "W": 186.2132,
    "Y": 163.1760,
    "V": 99.1326,
}
WATER_MASS = 18.01528

# Standard pKa values for charge titration and pI estimation
PKA_VALUES = {
    "N_TERM": 9.69,
    "C_TERM": 2.34,
    "C": 8.33,
    "D": 3.86,
    "E": 4.25,
    "H": 6.00,
    "K": 10.53,
    "R": 12.48,
    "Y": 10.07,
}

VALID_AMINO_ACIDS = set(AA_RESIDUE_MASSES.keys())


class ProteaseEngine:
    """
    High-fidelity peptide cleavage analyzer modeling physiological and diagnostic proteases.
    """

    def __init__(self) -> None:
        self._proteases: Dict[str, ProteaseSpec] = self._init_protease_catalog()

    @property
    def available_proteases(self) -> List[str]:
        """Return list of all registered protease identifiers."""
        return sorted(list(self._proteases.keys()))

    def get_protease_spec(self, name: str) -> Optional[ProteaseSpec]:
        """Retrieve metadata and specificity for a specific protease."""
        return self._proteases.get(name.lower())

    def _init_protease_catalog(self) -> Dict[str, ProteaseSpec]:
        """Initialize canonical protease specificity database."""
        catalog = {
            "trypsin": ProteaseSpec(
                name="Trypsin",
                family="Serine",
                description="Cleaves C-terminal to Lys or Arg, blocked by P1' Proline",
                biological_context="digestive / serum",
                optimum_ph=8.0,
                p1_allowed=["K", "R"],
                p1_prime_disallowed=["P"],
            ),
            "chymotrypsin_high": ProteaseSpec(
                name="Chymotrypsin (High specificity)",
                family="Serine",
                description="Cleaves C-terminal to Phe, Trp, Tyr, blocked by P1' Proline",
                biological_context="digestive",
                optimum_ph=7.8,
                p1_allowed=["F", "W", "Y"],
                p1_prime_disallowed=["P"],
            ),
            "chymotrypsin_low": ProteaseSpec(
                name="Chymotrypsin (Low specificity)",
                family="Serine",
                description="Cleaves C-terminal to Phe, Trp, Tyr, Leu, Met, blocked by P1' Proline",
                biological_context="digestive",
                optimum_ph=7.8,
                p1_allowed=["F", "W", "Y", "L", "M"],
                p1_prime_disallowed=["P"],
            ),
            "pepsin_ph13": ProteaseSpec(
                name="Pepsin (pH 1.3)",
                family="Aspartic",
                description="Cleaves C-terminal to Phe, Leu, Tyr, Trp at gastric pH",
                biological_context="digestive",
                optimum_ph=1.3,
                p1_allowed=["F", "L", "W", "Y"],
                p1_prime_disallowed=["P"],
            ),
            "pepsin_ph2": ProteaseSpec(
                name="Pepsin (pH > 2)",
                family="Aspartic",
                description="Broad hydrophobic cleavage at P1 (Phe, Leu, Met, Trp, Tyr, Glu)",
                biological_context="digestive",
                optimum_ph=2.0,
                p1_allowed=["F", "L", "M", "W", "Y", "E"],
                p1_prime_disallowed=["P"],
            ),
            "elastase": ProteaseSpec(
                name="Elastase",
                family="Serine",
                description="Cleaves C-terminal to small aliphatic residues (Ala, Val, Gly, Ser)",
                biological_context="digestive / neutrophil",
                optimum_ph=8.5,
                p1_allowed=["A", "V", "G", "S", "I", "L"],
                p1_prime_disallowed=["P"],
            ),
            "proteinase_k": ProteaseSpec(
                name="Proteinase K",
                family="Serine",
                description="Broad cleavage C-terminal to aliphatic and aromatic amino acids",
                biological_context="microbial / robust",
                optimum_ph=7.5,
                p1_allowed=["A", "E", "F", "I", "L", "T", "V", "W", "Y"],
                p1_prime_disallowed=["P"],
            ),
            "thermolysin": ProteaseSpec(
                name="Thermolysin",
                family="Metalloprotease",
                description="Cleaves N-terminal to bulky hydrophobic residues (Ile, Leu, Val, Phe, Met, Ala)",
                biological_context="microbial / experimental",
                optimum_ph=7.5,
                p1_prime_allowed=["I", "L", "V", "F", "M", "A"],
                p1_disallowed=["P"],
            ),
            "endoproteinase_gluc": ProteaseSpec(
                name="Endoproteinase GluC (V8)",
                family="Serine",
                description="Cleaves C-terminal to Glu (and Asp in phosphate buffer)",
                biological_context="bacterial",
                optimum_ph=7.8,
                p1_allowed=["E"],
                p1_prime_disallowed=["P"],
            ),
            "endoproteinase_lysc": ProteaseSpec(
                name="Endoproteinase LysC",
                family="Serine",
                description="Strictly cleaves C-terminal to Lys",
                biological_context="bacterial",
                optimum_ph=8.5,
                p1_allowed=["K"],
                p1_prime_disallowed=["P"],
            ),
            "endoproteinase_aspn": ProteaseSpec(
                name="Endoproteinase AspN",
                family="Metalloprotease",
                description="Cleaves N-terminal to Asp and Glu residues",
                biological_context="bacterial",
                optimum_ph=8.0,
                p1_prime_allowed=["D", "E"],
                p1_disallowed=["P"],
            ),
            "thrombin": ProteaseSpec(
                name="Thrombin",
                family="Serine",
                description="Cleaves Arg bonds with high preference for Gly at P1'",
                biological_context="serum / coagulation",
                optimum_ph=7.4,
                p1_allowed=["R"],
                p1_prime_allowed=["G"],
                p1_prime_disallowed=["P", "D", "E"],
            ),
            "factor_xa": ProteaseSpec(
                name="Factor Xa",
                family="Serine",
                description="Recognizes Ile-Glu/Asp-Gly-Arg motif and cleaves after Arg",
                biological_context="serum / coagulation",
                optimum_ph=7.4,
                p1_allowed=["R"],
                p1_prime_disallowed=["P"],
            ),
            "dpp4": ProteaseSpec(
                name="Dipeptidyl Peptidase IV (DPP-IV)",
                family="Dipeptidyl",
                description="Cleaves N-terminal dipeptide if penultimate residue is Pro or Ala",
                biological_context="serum / plasma",
                optimum_ph=7.4,
                n_terminal_only=True,
            ),
            "neprilysin": ProteaseSpec(
                name="Neprilysin (NEP)",
                family="Metalloprotease",
                description="Cleaves N-terminal to hydrophobic residues (Phe, Leu, Ile)",
                biological_context="kidney / plasma membrane",
                optimum_ph=7.0,
                p1_prime_allowed=["F", "L", "I"],
                p1_disallowed=["P"],
            ),
            "cathepsin_b": ProteaseSpec(
                name="Cathepsin B",
                family="Cysteine",
                description="Lysosomal cysteine protease, cleaves after Arg-Arg or Phe-Arg",
                biological_context="lysosomal",
                optimum_ph=5.0,
                p1_allowed=["R"],
                p1_prime_disallowed=["P"],
            ),
            "cathepsin_d": ProteaseSpec(
                name="Cathepsin D",
                family="Aspartic",
                description="Lysosomal aspartic protease, cleaves hydrophobic pairs (Phe, Leu)",
                biological_context="lysosomal",
                optimum_ph=4.5,
                p1_allowed=["F", "L", "Y"],
                p1_prime_allowed=["F", "L", "Y", "V", "I"],
                p1_prime_disallowed=["P"],
            ),
        }
        return catalog

    def validate_sequence(self, sequence: str) -> str:
        """Sanitize, uppercase, and validate amino acid sequence."""
        clean_seq = re.sub(r"\s+", "", sequence).upper()
        if not clean_seq:
            raise ValueError("Peptide sequence cannot be empty.")
        invalid_chars = set(clean_seq) - VALID_AMINO_ACIDS
        if invalid_chars:
            raise ValueError(
                f"Invalid non-canonical amino acid(s) detected: {sorted(list(invalid_chars))}. "
                f"Supported 20 standard IUPAC amino acids only."
            )
        return clean_seq

    def calculate_molecular_weight(self, sequence: str) -> float:
        """Compute average molecular mass in Daltons including terminal H and OH."""
        clean_seq = self.validate_sequence(sequence)
        total = sum(AA_RESIDUE_MASSES[aa] for aa in clean_seq) + WATER_MASS
        return total

    def calculate_isoelectric_point(self, sequence: str) -> float:
        """Estimate theoretical isoelectric point (pI) using bisection titration."""
        clean_seq = self.validate_sequence(sequence)

        # Count ionizable residues
        counts = {aa: clean_seq.count(aa) for aa in ["C", "D", "E", "H", "K", "R", "Y"]}

        def net_charge(ph: float) -> float:
            # Positive charges: N-terminus, Lys, Arg, His
            c_n_term = 1.0 / (1.0 + 10.0 ** (ph - PKA_VALUES["N_TERM"]))
            c_k = counts["K"] / (1.0 + 10.0 ** (ph - PKA_VALUES["K"]))
            c_r = counts["R"] / (1.0 + 10.0 ** (ph - PKA_VALUES["R"]))
            c_h = counts["H"] / (1.0 + 10.0 ** (ph - PKA_VALUES["H"]))
            pos = c_n_term + c_k + c_r + c_h

            # Negative charges: C-terminus, Asp, Glu, Cys, Tyr
            c_c_term = 1.0 / (1.0 + 10.0 ** (PKA_VALUES["C_TERM"] - ph))
            c_d = counts["D"] / (1.0 + 10.0 ** (PKA_VALUES["D"] - ph))
            c_e = counts["E"] / (1.0 + 10.0 ** (PKA_VALUES["E"] - ph))
            c_c = counts["C"] / (1.0 + 10.0 ** (PKA_VALUES["C"] - ph))
            c_y = counts["Y"] / (1.0 + 10.0 ** (PKA_VALUES["Y"] - ph))
            neg = c_c_term + c_d + c_e + c_c + c_y

            return pos - neg

        low = 0.0
        high = 14.0
        for _ in range(50):
            mid = (low + high) / 2.0
            charge = net_charge(mid)
            if charge > 0:
                low = mid
            else:
                high = mid
        return (low + high) / 2.0

    def _extract_cleavage_window(self, sequence: str, site_idx: int, span: int = 4) -> str:
        """Generate Schechter & Berger representation (e.g. P4P3P2P1|P1'P2'P3'P4')."""
        # site_idx is 1-based (bond between site_idx and site_idx + 1)
        start_left = max(0, site_idx - span)
        left = sequence[start_left:site_idx]

        end_right = min(len(sequence), site_idx + span)
        right = sequence[site_idx:end_right]

        return f"{left}|{right}"

    def analyze(
        self,
        sequence: str,
        peptide_id: str = "peptide_candidate",
        proteases: Optional[List[str]] = None,
    ) -> CleavageResult:
        """
        Perform comprehensive cleavage profiling of the peptide.
        """
        clean_seq = self.validate_sequence(sequence)
        seq_len = len(clean_seq)
        mw = self.calculate_molecular_weight(clean_seq)

        if proteases is None:
            active_protease_keys = self.available_proteases
        else:
            active_protease_keys = [p.lower() for p in proteases]
            # validate requested
            for pk in active_protease_keys:
                if pk not in self._proteases:
                    raise ValueError(f"Unknown protease: {pk}. Available: {self.available_proteases}")

        detected_sites: List[CleavageSite] = []
        cleavage_map: Dict[int, List[str]] = {i: [] for i in range(1, seq_len)}

        for pkey in active_protease_keys:
            spec = self._proteases[pkey]

            # Special cases:
            # 1. DPP-4: Cleaves dipeptide from N-terminus if residue 2 is Pro or Ala
            if pkey == "dpp4":
                if seq_len >= 3:
                    p1_res = clean_seq[1]  # 2nd residue (0-indexed 1)
                    if p1_res in ["P", "A"]:
                        site_idx = 2
                        w = self._extract_cleavage_window(clean_seq, site_idx)
                        site = CleavageSite(
                            site_index=site_idx,
                            protease_name=spec.name,
                            p1_residue=p1_res,
                            p1_prime_residue=clean_seq[2],
                            cleavage_window=w,
                            specificity_class="canonical",
                            cleavage_efficiency=1.0,
                        )
                        detected_sites.append(site)
                        cleavage_map[site_idx].append(spec.name)
                continue

            # 2. Factor Xa: Ile-Glu/Asp-Gly-Arg |
            if pkey == "factor_xa":
                for i in range(4, seq_len):  # need 4 residues before cleavage
                    p4 = clean_seq[i - 4]
                    p3 = clean_seq[i - 3]
                    p2 = clean_seq[i - 2]
                    p1 = clean_seq[i - 1]
                    p1_prime = clean_seq[i]
                    if p4 == "I" and p3 in ["E", "D"] and p2 == "G" and p1 == "R" and p1_prime != "P":
                        site_idx = i
                        w = self._extract_cleavage_window(clean_seq, site_idx)
                        site = CleavageSite(
                            site_index=site_idx,
                            protease_name=spec.name,
                            p1_residue=p1,
                            p1_prime_residue=p1_prime,
                            cleavage_window=w,
                            specificity_class="canonical",
                            cleavage_efficiency=1.0,
                        )
                        detected_sites.append(site)
                        cleavage_map[site_idx].append(spec.name)
                continue

            # Standard scanning across all internal peptide bonds (1 to seq_len - 1)
            for i in range(1, seq_len):
                p1_res = clean_seq[i - 1]
                p1_prime_res = clean_seq[i]

                # P1' check
                if spec.p1_prime_disallowed and p1_prime_res in spec.p1_prime_disallowed:
                    continue
                if spec.p1_prime_allowed and p1_prime_res not in spec.p1_prime_allowed:
                    continue

                # P1 check
                if spec.p1_disallowed and p1_res in spec.p1_disallowed:
                    continue
                if spec.p1_allowed and p1_res not in spec.p1_allowed:
                    continue

                # Specific enzyme subsite filters
                if pkey == "thrombin":
                    # Thrombin requires Arg at P1 and Gly at P1' with high selectivity
                    if p1_res != "R" or p1_prime_res != "G":
                        continue

                if pkey == "cathepsin_b":
                    # Cathepsin B prefers basic/hydrophobic at P2: Arg-Arg or Phe-Arg
                    if i >= 2:
                        p2_res = clean_seq[i - 2]
                        if p2_res not in ["R", "F", "V", "L"]:
                            continue
                    else:
                        continue

                # Cleavage confirmed
                site_idx = i
                w = self._extract_cleavage_window(clean_seq, site_idx)
                site = CleavageSite(
                    site_index=site_idx,
                    protease_name=spec.name,
                    p1_residue=p1_res,
                    p1_prime_residue=p1_prime_res,
                    cleavage_window=w,
                    specificity_class="canonical",
                    cleavage_efficiency=1.0,
                )
                detected_sites.append(site)
                cleavage_map[site_idx].append(spec.name)

        # Generate digestion fragments from unique cut positions
        cut_positions = sorted(list({s.site_index for s in detected_sites}))
        fragments: List[PeptideFragment] = []
        prev_cut = 0

        for frag_idx, cut in enumerate(cut_positions + [seq_len], start=1):
            frag_seq = clean_seq[prev_cut:cut]
            if frag_seq:
                frag_mw = sum(AA_RESIDUE_MASSES[aa] for aa in frag_seq) + WATER_MASS
                frag_pi = self.calculate_isoelectric_point(frag_seq)
                fragments.append(
                    PeptideFragment(
                        fragment_index=frag_idx,
                        start_pos=prev_cut + 1,
                        end_pos=cut,
                        sequence=frag_seq,
                        length=len(frag_seq),
                        molecular_weight=frag_mw,
                        isoelectric_point=frag_pi,
                    )
                )
            prev_cut = cut

        # Filter out empty map entries
        active_cleavage_map = {k: v for k, v in cleavage_map.items() if v}

        # Calculate vulnerability index & resistance score
        total_cuts = len(detected_sites)
        possible_bonds = max(1, seq_len - 1)
        num_proteases = max(1, len(active_protease_keys))

        # Normalized degradation load per bond and per protease
        degradation_load = total_cuts / (possible_bonds * num_proteases)
        pvi = min(100.0, degradation_load * 100.0 * 5.0)  # scaled
        resistance_score = max(0.0, 100.0 - pvi)

        # Identify hot spots (residues with highest cleavage multiplicity)
        site_counts = {k: len(v) for k, v in active_cleavage_map.items()}
        sorted_hotspots = sorted(site_counts.keys(), key=lambda k: site_counts[k], reverse=True)
        hotspot_residues = [pos for pos in sorted_hotspots if site_counts[pos] >= 2]

        return CleavageResult(
            peptide_id=peptide_id,
            sequence=clean_seq,
            length=seq_len,
            molecular_weight=mw,
            proteases_analyzed=[self._proteases[k].name for k in active_protease_keys],
            cleavage_sites=detected_sites,
            fragments=fragments,
            cleavage_count=total_cuts,
            proteolytic_vulnerability_index=pvi,
            proteolytic_resistance_score=resistance_score,
            hotspot_residues=hotspot_residues,
            cleavage_map=active_cleavage_map,
        )
