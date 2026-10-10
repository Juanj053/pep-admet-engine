rule plasma_halflife:
    input:
        cleavage_json="results/{sample}/cleavage_profile.json",
    output:
        halflife_json=report(
            "results/{sample}/plasma_halflife.json",
            caption="../report/plasma_halflife.rst",
            category="Pharmacokinetics & Distribution",
            subcategory="Systemic Half-Life",
        ),
    log:
        "logs/{sample}/plasma_halflife.log",
    conda:
        "../envs/plasma_halflife.yaml"
    params:
        workflow_basedir=workflow.basedir,
        sequence=get_peptide_sequence,
    script:
        "../scripts/run_plasma_halflife.py"
