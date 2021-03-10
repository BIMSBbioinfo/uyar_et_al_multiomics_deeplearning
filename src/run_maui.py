import argparse
import os
import pandas as pd
import maui
from arcaslib import squishy_maui

parser = argparse.ArgumentParser(description='Run MAUI to get latent factors from the given assay matrix.')

parser.add_argument("assay_matrix", help = "Path to CSV format feature x sample matrix")
parser.add_argument("outFile", help = "Path to output file for latent factors")
parser.add_argument("Nlatent", help = "Number of latent factors to search for")
parser.add_argument("Nepochs", help = "Number of epochs to run the training algorithm")
parser.add_argument("Nthreads", help = "Number of threads used to run maui")
args = parser.parse_args()

#read assay matrix
assay = pd.read_csv(args.assay_matrix) 
print(assay.head())
#start session
from keras import backend as K
K.set_session(K.tf.Session(config=K.tf.ConfigProto(intra_op_parallelism_threads=int(args.Nthreads), inter_op_parallelism_threads=int(args.Nthreads))))
maui_model = maui.Maui(n_hidden=[1000], n_latent=int(args.Nlatent), epochs=int(args.Nepochs))
z = maui_model.fit_transform({'assay': assay})
# save LF to file
z.to_csv(args.outFile)

# get feature contributions to each LF
nwp = maui_model.get_neural_weight_product()
nwp.to_csv(''.join([os.path.splitext(args.outFile)[0], '.feature_weights.csv']))
