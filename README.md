


# Downloading TCGA

## Omics Data

### GEX + MUT + METH: download.tcga.R
This script uses TCGAbiolinks package to download gene expression (FPKM), methylation, and mutation data
files along with patient-related clinical meta-data. 

The data is first downloaded into the "GDCdata" sub-folder in the current working directory. 
From there the data is further processed and prepared, which is written under "GDCdata_prepared".

usage: 
<path to Rscript> ./src/download.tcga.R <target folder to save GDCdata> <target folder to save prepared data>
example: 
/opt/R/4.0/bin/Rscript ./src/download.tcga.R `pwd` `pwd`/GDCdata_prepared

### CNVs: download.tcga.firehore.R

This script uses RTCGAToolbox and TCGAbiolinks packages to download Somatic Copy Number Alteration (GISTIC scores) 
from BROAD Firehose data. 

The data is first downloaded into the "/GDCdata/FirehoseData" sub-folder in the current working directory.
From there the data is further processed and prepared, which is written under "GDCdata_prepared".

usage:
/opt/R/4.0/bin/Rscript ./src/download.tcga.firehose.R `pwd`/GDCdata/FirehoseData `pwd`/GDCdata_prepared 

### MSI: download.tcga.MSI.R

This script uses TCGAbiolinks data to download Microsatellite Instability (MSI) status annotations for 
TCGA projects (pan-gastrointestinal cancers). 

The data is first downloaded into the "/GDCdata/FirehoseData" sub-folder in the current working directory.
From there the data is further processed and prepared, which is written under "GDCdata_prepared".

usage: 
/opt/R/4.0/bin/Rscript ./src/download.tcga.MSI.R `pwd` `pwd`/GDCdata_prepared

### Clinical data

Pancancer survival and further clinical data from [Liu J. et al, Cell, 2018](https://www.sciencedirect.com/science/article/pii/S0092867418302290?via%3Dihub)
can be imported using the `RDS` object under `./data/TCGA.surv.RDS`.

Pancancer clinical data downloaded from TCGA can be imported using the `RDS` object under `./data/TCGA.clin.RDS`. 

# Analysis

Here we describe how to prepare the omics datasets and how to run multi-omics integration tools: MAUI, MOFA, MCIA, and PCA. 

`setup_experiments.R` script prepares the inputs and `snakefile.py` script runs a snakemake pipeline on all the inputs. 

Both scripts require a `settings.yaml` file which includes all necessary input and tool configurations. 

See `./settings.yaml` file that was used as input for both scripts.

## Preparing omics datasets 

Here we prepare the input files for multi-omics integration tools. The input files for 
these tools will be prepared under `./assays` folder. 

usage: 
/opt/R/4.0/bin/Rscript ./src/setup_experiments.R --settings ./settings.yaml

## Running multi-omics integration tools 

Here we run a snakemake pipeline which invokes commands to run MAUI, MOFA, MCIA, and PCA to 
do multi-omics integration. The outputs including the learned latent factors along with feature importance 
values are written under `./output`. 

usage:
time snakemake -p -s ./src/snakefile.py -j 10 --configfile ./settings.yaml --keep-going

# Manuscript Figures

`compile_figures.sh` script is used to invoke `manuscript_figures.Rmd` file to make the manuscript figures in this study. 

usage: 

manuscript figures:
# for releasing
nohup time -v bash ./src/compile_figures.sh /opt/R/4.0/bin/Rscript ./src/manuscript_figures.Rmd `pwd`/results/hallmarks_xcell/settings.yaml /data/bimsbstatic/public/akalin/buyar/manuscript_figures_arcas /data/bimsbstatic/public/akalin/buyar/manuscript_figures_arcas > figures.log
