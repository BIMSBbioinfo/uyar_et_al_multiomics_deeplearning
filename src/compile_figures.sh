RSCRIPT=$1
RMD=$2 # path the Rmd file
settings_file=$3 # settings file used to run the pipeline
figures_folder=$4  # where to export pdf figures
workdir=$5 # where to save docx file

${RSCRIPT} -e "library(rmarkdown); rmarkdown::render(\"${RMD}\", output_dir = \"${workdir}\", params = list(settings_file = \"${settings_file}\", figures_folder = \"${figures_folder}\", workdir = \"${workdir}\"))"
