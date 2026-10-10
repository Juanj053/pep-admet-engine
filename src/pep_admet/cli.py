"""
PEP-ADMET Engine Command Line Interface.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from pep_admet.core.halflife_predictor import PlasmaHalfLifePredictor
from pep_admet.core.permeability_analyzer import MembranePermeabilityAnalyzer
from pep_admet.core.protease_engine import ProteaseEngine
from pep_admet.core.report_builder import ADMETReportBuilder
from pep_admet.core.stabilization_advisor import StabilizationAdvisor


def main() -> int:
    parser = argparse.ArgumentParser(
        description="PEP-ADMET Engine: In-Silico ADMET Profiling and Rational Stabilization for Therapeutic Peptides."
    )
    parser.add_argument("--sequence", "-s", type=str, required=True, help="Peptide amino acid sequence")
    parser.add_argument("--name", "-n", type=str, default="peptide_candidate", help="Peptide candidate identifier")
    parser.add_argument("--output", "-o", type=str, default=None, help="Path to output JSON profile")
    args = parser.parse_args()

    protease_eng = ProteaseEngine()
    halflife_pred = PlasmaHalfLifePredictor()
    perm_analyzer = MembranePermeabilityAnalyzer()
    stab_advisor = StabilizationAdvisor()
    builder = ADMETReportBuilder()

    cleavage = protease_eng.analyze(args.sequence, peptide_id=args.name)
    halflife = halflife_pred.predict(args.sequence, peptide_id=args.name, proteolytic_vulnerability_index=cleavage.proteolytic_vulnerability_index)
    perm = perm_analyzer.analyze(args.sequence, peptide_id=args.name)
    plan = stab_advisor.formulate_plan(cleavage, halflife, perm)
    profile = builder.build_profile(cleavage, halflife, perm, plan)

    out_data = profile.to_dict()
    if args.output:
        builder.export_json(profile, args.output)
        print(f"Profile written to {args.output}")
    else:
        print(json.dumps(out_data, indent=2))

    return 0


if __name__ == "__main__":
    sys.exit(main())
