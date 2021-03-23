Run pan-cancer, pan-organ experiments and compare with other tools

1. Set up the experiments

Configure the experiment using the settings.yaml
> /opt/R-4.0.2/lib64/R/bin/Rscript ./src/setup_experiments.R --settings ./settings.yaml

2. Run pipeline

> time snakemake -p -s ./src/snakefile.py -j 10 --configfile ./settings2.yaml --keep-going

./results folder contains sub-folders, where each sub-folder is for a different kind of analysis 
Each sub-folder should contain an "assay" folder, a "settings.yaml" file for that experiment, and an "output" folder for the pipeline's output.


