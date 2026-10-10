rule reporting:
    input:
        cleavage_jsons=expand("results/{sample}/cleavage_profile.json", sample=config["peptides"]),
        halflife_jsons=expand("results/{sample}/plasma_halflife.json", sample=config["peptides"]),
        permeability_jsons=expand("results/{sample}/membrane_permeability.json", sample=config["peptides"]),
        stabilization_jsons=expand("results/{sample}/stabilization_plan.json", sample=config["peptides"]),
    output:
        per_sample_admet=expand("results/{sample}/admet_profile.json", sample=config["peptides"]),
        per_sample_plots=expand("results/{sample}/admet_analytics.png", sample=config["peptides"]),
        summary_csv=report(
            "results/admet_summary.csv",
            caption="../report/reporting.rst",
            category="Synthesis & Summary",
            subcategory="Comparative Metrics Table",
        ),
        global_report_json=report(
            "results/global_admet_report.json",
            caption="../report/reporting.rst",
            category="Synthesis & Summary",
            subcategory="Global JSON Export",
        ),
    log:
        "logs/reporting.log",
    conda:
        "../envs/reporting.yaml"
    params:
        workflow_basedir=workflow.basedir,
        dpi=config.get("dpi", 300),
        generate_plots=config.get("generate_plots", True),
    script:
        "../scripts/run_reporting.py"
