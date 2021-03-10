Run pan-cancer, pan-organ experiments and compare with other tools

1. Set up the experiments
Configure the experiment using the settings.yaml
> /opt/R-4.0.2/lib64/R/bin/Rscript ./src/setup_experiments.R --settings ./settings.yaml

2. Run pipeline
# modify the output path in snakefile.py if needed
# using 20 cores
nohup time snakemake -p -s ./src/snakefile.py -j 10 --keep-going
