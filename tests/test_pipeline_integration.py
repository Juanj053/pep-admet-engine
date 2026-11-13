"""
Integration tests for PEP-ADMET complete workflow.
"""

import json
import sys

from pep_admet import cli
from pep_admet.core import (
    ADMETReportBuilder,
    MembranePermeabilityAnalyzer,
    PlasmaHalfLifePredictor,
    ProteaseEngine,
    StabilizationAdvisor,
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


def test_compare_analogs_and_cli(tmp_path, monkeypatch):
    seq = "HAEGTFTSDVSSYLEGQAAKEFIAWLVKGR"
    protease_eng = ProteaseEngine()
    halflife_pred = PlasmaHalfLifePredictor()
    perm_analyzer = MembranePermeabilityAnalyzer()
    stab_advisor = StabilizationAdvisor()
    builder = ADMETReportBuilder()

    cleavage = protease_eng.analyze(seq, peptide_id="glp1_native")
    halflife = halflife_pred.predict(
        seq,
        peptide_id="glp1_native",
        proteolytic_vulnerability_index=cleavage.proteolytic_vulnerability_index,
    )
    perm = perm_analyzer.analyze(seq, peptide_id="glp1_native")
    plan = stab_advisor.formulate_plan(cleavage, halflife, perm)
    profile = builder.build_profile(cleavage, halflife, perm, plan)

    comparisons = builder.compare_analogs(profile, protease_eng, halflife_pred, perm_analyzer)
    assert len(comparisons) >= 3
    required_keys = {
        "variant_id",
        "sequence",
        "proteolytic_vulnerability_index",
        "estimated_t_half_minutes",
        "half_life_fold_vs_wt",
        "caco2_log_papp_cm_s",
        "overall_admet_score",
    }
    for row in comparisons:
        assert required_keys.issubset(row.keys())
        assert float(row["estimated_t_half_minutes"]) > 0.0

    out_file = tmp_path / "glp1_analogs.json"
    monkeypatch.setattr(
        sys,
        "argv",
        ["pep-admet", "--sequence", seq, "--name", "glp1", "--compare-analogs", "--output", str(out_file)],
    )
    assert cli.main() == 0

    payload = json.loads(out_file.read_text(encoding="utf-8"))
    assert "analog_comparison" in payload
    assert len(payload["analog_comparison"]) == len(comparisons)
