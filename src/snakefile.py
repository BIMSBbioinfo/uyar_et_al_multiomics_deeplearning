
import os

ASSAYDIR = config['assay_output']['folder'] 
SRCDIR = '/data/local/buyar/arcas/pancancer_multiomics_manuscript/src' 
ASSAYS = set([os.path.splitext(f)[0] for f in os.listdir(ASSAYDIR) if re.match(r'.*\.csv$', f)])
RSCRIPT = '/opt/R-4.0.2/lib64/R/bin/Rscript'
OUTDIR = config['pipeline_output']
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
        script = os.path.join(SRCDIR, 'run_MCIA.R'), 
        nfactors = config['mcia']['nfactors']
    shell:
        "{RSCRIPT} {params.script} {input} {output[0]} {params.nfactors} > {log} 2>&1"

rule run_maui:
    input:
        os.path.join(ASSAYDIR, "{assay}.csv")
    output:
        os.path.join(OUTDIR, 'maui_factors', "{assay}.maui_factors.csv"),
        os.path.join(OUTDIR, 'maui_factors', "{assay}.maui_factors.feature_weights.csv")
    log: os.path.join(LOG_DIR, "run_maui.{assay}.log")
    params:
        script = os.path.join(SRCDIR, 'run_maui.sh'),
        nfactors = config['maui']['nfactors'],
        epochs = config['maui']['epochs'],
        threads = config['maui']['threads']
    shell:
        "bash {params.script} {input} {output[0]} {params.nfactors} {params.epochs} {params.threads} > {log} 2>&1"

rule run_mofa:
    input:
        os.path.join(ASSAYDIR, "{assay}.RDS")
    output:
        os.path.join(OUTDIR, 'mofa_factors', "{assay}.mofa_factors.csv"),
        os.path.join(OUTDIR, 'mofa_factors', "{assay}.mofa_factors.feature_weights.csv")
    log: os.path.join(LOG_DIR, "run_mofa.{assay}.log")
    params:
        script = os.path.join(SRCDIR, 'run_mofa.R'),
        nfactors = config['mofa']['nfactors']
    shell:
        "export OMP_NUM_THREADS=2; {RSCRIPT} {params.script} {input} {output[0]} {params.nfactors} > {log} 2>&1"

rule run_pca:
    input:
        os.path.join(ASSAYDIR, "{assay}.csv")
    output:
        os.path.join(OUTDIR, 'pca_factors', "{assay}.pca_factors.csv"),
        os.path.join(OUTDIR, 'pca_factors', "{assay}.pca_factors.feature_weights.csv")
    log: os.path.join(LOG_DIR, "run_pca.{assay}.log")
    params:
        script = os.path.join(SRCDIR, 'run_PCA.R'),
        nfactors = config['pca']['nfactors']
    shell:
        "{RSCRIPT} {params.script} {input} {output[0]} {params.nfactors} > {log} 2>&1"

