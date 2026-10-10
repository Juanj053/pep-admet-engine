"""
Unit tests for ProteaseEngine.
"""

import pytest
import sys
from pathlib import Path

# Add src to path
src_dir = Path(__file__).resolve().parent.parent / 'src'
if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

from pep_admet.core.protease_engine import ProteaseEngine
from pep_admet.models.cleavage import CleavageResult, CleavageSite


@pytest.fixture
def engine():
    return ProteaseEngine()


def test_protease_catalog_initialization(engine):
    proteases = engine.available_proteases
    assert "trypsin" in proteases
    assert "chymotrypsin_high" in proteases
    assert "pepsin_ph13" in proteases
    assert "dpp4" in proteases
    assert "thrombin" in proteases
    assert len(proteases) >= 15


def test_sequence_validation(engine):
    # Valid uppercase sequence
    assert engine.validate_sequence("ACDEFGHIKLMNPQRSTVWY") == "ACDEFGHIKLMNPQRSTVWY"
    # Whitespace and lowercase
    assert engine.validate_sequence(" acd  efg ") == "ACDEFG"

    # Non-canonical amino acid
    with pytest.raises(ValueError, match="Invalid non-canonical amino acid"):
        engine.validate_sequence("ACDEFGHIKX")

    # Empty sequence
    with pytest.raises(ValueError, match="cannot be empty"):
        engine.validate_sequence("   ")


def test_molecular_weight_and_pi(engine):
    # Dipeptide Ala-Ala: MW = 71.0788 * 2 + 18.01528 = 160.17
    mw = engine.calculate_molecular_weight("AA")
    assert abs(mw - 160.17) < 0.1

    pi = engine.calculate_isoelectric_point("AA")
    assert 5.0 < pi < 7.0


def test_trypsin_cleavage_and_proline_block(engine):
    # Sequence with K and R: "AGKAGR" -> cuts after K(3) and R(6 is C-term so no internal cut)
    res = engine.analyze("AGKAGR", proteases=["trypsin"])
    assert res.cleavage_count == 1
    assert res.cleavage_sites[0].site_index == 3
    assert res.cleavage_sites[0].p1_residue == "K"
    assert res.cleavage_sites[0].p1_prime_residue == "A"

    # Proline block: "AGKPAGR" -> K is followed by P, so no cleavage!
    res_blocked = engine.analyze("AGKPAGR", proteases=["trypsin"])
    # Cuts: bond 4 is P-A, bond 6 is G-R. K is at pos 3, followed by P at pos 4. Blocked!
    sites_at_k = [s for s in res_blocked.cleavage_sites if s.site_index == 3]
    assert len(sites_at_k) == 0


def test_dpp4_n_terminal_dipeptide_cleavage(engine):
    # GLP-1 starts with HAE... (His-Ala-Glu): Ala at pos 2 triggers DPP-4
    res = engine.analyze("HAEGTFTSD", proteases=["dpp4"])
    assert res.cleavage_count == 1
    assert res.cleavage_sites[0].site_index == 2
    assert res.cleavage_sites[0].p1_residue == "A"

    # Substance P starts with RPK... (Arg-Pro-Lys): Pro at pos 2 triggers DPP-4
    res_sp = engine.analyze("RPKPQQFFGLM", proteases=["dpp4"])
    assert res_sp.cleavage_count == 1
    assert res_sp.cleavage_sites[0].site_index == 2
    assert res_sp.cleavage_sites[0].p1_residue == "P"

    # Non-DPP4 substrate: starts with HGE... (His-Gly-Glu) -> no cleavage
    res_none = engine.analyze("HGEGTFTSD", proteases=["dpp4"])
    assert res_none.cleavage_count == 0


def test_fragment_reconstruction(engine):
    # Sequence: "AGKFL"
    # Trypsin cleaves after K(3) -> fragments: [AGK, FL]
    res = engine.analyze("AGKFL", proteases=["trypsin"])
    assert len(res.fragments) == 2
    assert res.fragments[0].sequence == "AGK"
    assert res.fragments[0].start_pos == 1
    assert res.fragments[0].end_pos == 3
    assert res.fragments[1].sequence == "FL"
    assert res.fragments[1].start_pos == 4
    assert res.fragments[1].end_pos == 5

    # Total length matches
    assert sum(f.length for f in res.fragments) == len("AGKFL")


def test_full_proteolytic_profile(engine):
    # Substance P: RPKPQQFFGLM
    res = engine.analyze("RPKPQQFFGLM", peptide_id="substance_p")
    assert res.length == 11
    assert res.cleavage_count > 0
    assert 0.0 <= res.proteolytic_vulnerability_index <= 100.0
    assert 0.0 <= res.proteolytic_resistance_score <= 100.0
    assert isinstance(res.hotspot_residues, list)
    assert isinstance(res.cleavage_map, dict)
    d = res.to_dict()
    assert d["peptide_id"] == "substance_p"
