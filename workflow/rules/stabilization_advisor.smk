rule stabilization_advisor:
    input:
        cleavage_json="results/{sample}/cleavage_profile.json",
        halflife_json="results/{sample}/plasma_halflife.json",
        permeability_json="results/{sample}/membrane_permeability.json",
    output:
        stabilization_json=report(
            "results/{sample}/stabilization_plan.json",
            caption="../report/stabilization_advisor.rst",
            category="Lead Optimization",
            subcategory="Peptidomimetic Engineering Plan",
        ),
    log:
        "logs/{sample}/stabilization_advisor.log",
    conda:
        "../envs/stabilization_advisor.yaml"
    params:
        workflow_basedir=workflow.basedir,
    script:
        "../scripts/run_stabilization_advisor.py"
