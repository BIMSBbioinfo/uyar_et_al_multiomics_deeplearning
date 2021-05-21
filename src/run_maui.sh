source activate maui
export PYTHONPATH=/home/buyar/.conda/envs/maui/lib/python3.6/site-packages
assays=$1 # assay matrix file path
outFile=$2 # path to output file
N_lfs=$3 #number of latent factors
N_epochs=$4 #nmuber of epochs
N_threads=$5 # number of threads
batchSize=$6 # size of batches to feed to the network
python /data/local/buyar/arcas/subtyping_paper/src/run_maui.py ${assays} ${outFile} ${N_lfs} ${N_epochs} ${N_threads} ${batchSize}
