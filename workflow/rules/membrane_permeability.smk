rule membrane_permeability:
    output:
        permeability_json=report(
            "results/{sample}/membrane_permeability.json",
            caption="../report/membrane_permeability.rst",
            category="Absorption & Permeability",
            subcategory="Physicochemical Descriptors",
        ),
    log:
        "logs/{sample}/membrane_permeability.log",
    conda:
        "../envs/membrane_permeability.yaml"
    params:
        workflow_basedir=workflow.basedir,
        sequence=get_peptide_sequence,
    script:
        "../scripts/run_membrane_permeability.py"
