library(rmarkdown)

args = commandArgs(trailingOnly = TRUE)

RMD=args[1] # path the Rmd file
settings_file=args[2] # settings file used to run the pipeline
immther_projdir=args[3] # path to folder with immunotherapy analysis
tmz_gbm_figure_data_folder=args[4] # path to folder with tmz-gbm figure data
figures_folder=args[5]  # where to export pdf figures
workdir=args[6] # where to save docx file

rmarkdown::render(input = RMD, output_dir = workdir, 
                  params = list('settings_file' = settings_file, 
                                'immther_projdir' = immther_projdir, 
                                'tmz_gbm_figure_data_folder' = tmz_gbm_figure_data_folder, 
                                'figures_folder' = figures_folder,
                                'workdir' = workdir))
