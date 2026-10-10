rule protease_cleavage:
    output:
        cleavage_json=report(
            "results/{sample}/cleavage_profile.json",
            caption="../report/protease_cleavage.rst",
            category="Proteolytic Metabolism",
            subcategory="Cleavage Profile",
        ),
        fragments_csv=report(
            "results/{sample}/digestion_fragments.csv",
            caption="../report/protease_cleavage.rst",
            category="Proteolytic Metabolism",
            subcategory="Fragments Table",
        ),
    log:
        "logs/{sample}/protease_cleavage.log",
    conda:
        "../envs/protease_cleavage.yaml"
    params:
        workflow_basedir=workflow.basedir,
        sequence=get_peptide_sequence,
        proteases=config.get("proteases", None),
    script:
        "../scripts/run_protease_cleavage.py"
