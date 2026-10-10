"""
Integration tests for PEP-ADMET complete workflow.
"""

import sys
from pathlib import Path
import pytest

src_dir = Path(__file__).resolve().parent.parent / 'src'
if str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

from pep_admet.core import (
    ProteaseEngine,
    PlasmaHalfLifePredictor,
    MembranePermeabilityAnalyzer,
    StabilizationAdvisor,
    ADMETReportBuilder,
)


def test_full_admet_pipeline_semaglutide_backbone():
    # GLP-1 native sequence: HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR
    seq = "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR"
    
    protease_eng = ProteaseEngine()
    halflife_pred = PlasmaHalfLifePredictor()
    perm_analyzer = MembranePermeabilityAnalyzer()
    stab_advisor = StabilizationAdvisor()
    builder = ADMETReportBuilder()

    cleavage = protease_eng.analyze(seq, peptide_id="glp1_native")
    assert cleavage.cleavage_count > 0
    assert cleavage.length == 30

    halflife = halflife_pred.predict(
        seq,
        peptide_id="glp1_native",
        proteolytic_vulnerability_index=cleavage.proteolytic_vulnerability_index,
    )
    assert halflife.estimated_t_half_minutes > 0.0

    perm = perm_analyzer.analyze(seq, peptide_id="glp1_native")
    assert perm.tpsa > 500  # large peptide
    assert perm.caco2_log_papp_cm_s < -5.0

    plan = stab_advisor.formulate_plan(cleavage, halflife, perm)
    assert len(plan.recommendations) >= 3
    # Verify DPP-4 recommendation is present
    has_dpp4_rec = any("DPP" in r.rationale for r in plan.recommendations)
    assert has_dpp4_rec is True

    profile = builder.build_profile(cleavage, halflife, perm, plan)
    assert 0.0 <= profile.overall_admet_score <= 100.0
    d = profile.to_dict()
    assert "overall_admet_score" in d
    assert "cleavage" in d
    assert "halflife" in d
    assert "permeability" in d
    assert "stabilization" in d
