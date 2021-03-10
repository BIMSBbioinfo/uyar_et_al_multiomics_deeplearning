
import os
ASSAYDIR = '/data/local/buyar/arcas/pancancer_multiomics_manuscript/assays' 
SRCDIR = '/data/local/buyar/arcas/pancancer_multiomics_manuscript/src' 
ASSAYS = set([os.path.splitext(f)[0] for f in os.listdir(ASSAYDIR) if re.match(r'.*\.csv$', f)])
RSCRIPT = '/opt/R-4.0.2/lib64/R/bin/Rscript'
OUTDIR = 'output'
LOG_DIR = os.path.join(OUTDIR, "logs")

rule all:
    input:
        expand(os.path.join(OUTDIR, 'maui_factors', "{assay}.maui_factors.csv"), assay = ASSAYS),
        expand(os.path.join(OUTDIR, 'maui_factors', "{assay}.maui_factors.feature_weights.csv"), assay = ASSAYS),
        expand(os.path.join(OUTDIR, 'mofa_factors', "{assay}.mofa_factors.csv"), assay = ASSAYS),
        expand(os.path.join(OUTDIR, 'mofa_factors', "{assay}.mofa_factors.feature_weights.csv"), assay = ASSAYS),
        expand(os.path.join(OUTDIR, 'pca_factors', "{assay}.pca_factors.csv"), assay = ASSAYS),
        expand(os.path.join(OUTDIR, 'pca_factors', "{assay}.pca_factors.feature_weights.csv"), assay = ASSAYS),
	expand(os.path.join(OUTDIR, 'mcia_factors', "{assay}.mcia_factors.csv"), assay = ASSAYS), 
	expand(os.path.join(OUTDIR, 'mcia_factors', "{assay}.mcia_factors.feature_weights.csv"), assay = ASSAYS)

rule run_mcia:
    input:
        os.path.join(ASSAYDIR, "{assay}.RDS")
    output:
        os.path.join(OUTDIR, 'mcia_factors', "{assay}.mcia_factors.csv"),
        os.path.join(OUTDIR, 'mcia_factors', "{assay}.mcia_factors.feature_weights.csv")
    log: os.path.join(LOG_DIR, "run_mcia.{assay}.log")
    params:
        script = os.path.join(SRCDIR, 'run_MCIA.R')
    shell:
        "{RSCRIPT} {params.script} {input} {output[0]} > {log} 2>&1"

rule run_maui:
    input:
        os.path.join(ASSAYDIR, "{assay}.csv")
    output:
        os.path.join(OUTDIR, 'maui_factors', "{assay}.maui_factors.csv"),
        os.path.join(OUTDIR, 'maui_factors', "{assay}.maui_factors.feature_weights.csv")
    log: os.path.join(LOG_DIR, "run_maui.{assay}.log")
    params:
        script = os.path.join(SRCDIR, 'run_maui.sh')
    shell:
        "bash {params.script} {input} {output[0]} 100 500 2 > {log} 2>&1"

rule run_mofa:
    input:
        os.path.join(ASSAYDIR, "{assay}.RDS")
    output:
        os.path.join(OUTDIR, 'mofa_factors', "{assay}.mofa_factors.csv"),
        os.path.join(OUTDIR, 'mofa_factors', "{assay}.mofa_factors.feature_weights.csv")
    log: os.path.join(LOG_DIR, "run_mofa.{assay}.log")
    params:
        script = os.path.join(SRCDIR, 'run_mofa.R')
    shell:
        "export OMP_NUM_THREADS=2; {RSCRIPT} {params.script} {input} {output[0]}  > {log} 2>&1"

rule run_pca:
    input:
        os.path.join(ASSAYDIR, "{assay}.csv")
    output:
        os.path.join(OUTDIR, 'pca_factors', "{assay}.pca_factors.csv"),
        os.path.join(OUTDIR, 'pca_factors', "{assay}.pca_factors.feature_weights.csv")
    log: os.path.join(LOG_DIR, "run_pca.{assay}.log")
    params:
        script = os.path.join(SRCDIR, 'run_PCA.R')
    shell:
        "{RSCRIPT} {params.script} {input} {output[0]} > {log} 2>&1"

