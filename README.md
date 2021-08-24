Run pan-cancer, pan-organ experiments and compare with other tools

1. Set up the experiments

Configure the experiment using the settings.yaml # set up a folder under results folder 
> mkdir results/my_experiment; cd results/my_experiment
> /opt/R/4.0/bin/Rscript ../../src/setup_experiments.R --settings ./settings.yaml

2. Run pipeline

> time snakemake -p -s ./src/snakefile.py -j 10 --configfile ./settings2.yaml --keep-going

./results folder contains sub-folders, where each sub-folder is for a different kind of analysis 
Each sub-folder should contain an "assay" folder, a "settings.yaml" file for that experiment, and an "output" folder for the pipeline's output.

3. compile figures
bash ./src/compile_figures.sh /opt/R-4.0.2/lib64/R/bin/Rscript `pwd`/src/manuscript_figures.Rmd `pwd`/<path to settings> <path to where to save the figures> <path to where to save docx>


manuscript figures:
nohup bash ./src/compile_figures.sh /opt/R/4.0/bin/Rscript ./src/manuscript_figures.Rmd `pwd`/results/hallmarks_many_LFs/settings.yaml /data/bimsbstatic/public/akalin/buyar/manuscript_figures_arcas/ /data/bimsbstatic/public/akalin/buyar/manuscript_figures_arcas > figures.log



