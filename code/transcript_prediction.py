#Python script for running Predict.py
#Part of TWAS_for_protein Pipeline
#Created by Henry Wittich
'''This python script takes genotype dosage data from a population in the TOPMed cohort, along with a GTEx
tissue expression model and runs PrediXcan to impute expression levels from individual genotypes, producing
the following output file:

tissue_type_predict.txt
tissue_type_summary.txt'''

import argparse
import sys
import os

cwd = os.getcwd() #current working directory

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Script to run Predict.py')
    parser.add_argument('--dose',
                        help='path to the dosage files',
                        required='True')
    parser.add_argument('--exp',
                        help='path to the expression model',
                        required='True')
    parser.add_argument('--samples',
                        help='path to the list of samples',
                        required='True')
    parser.add_argument('--out','-o',
                        help='path to output file',
                        default=cwd)
    return parser.parse_args(args)

#retrieve command line arguments
args=check_arg(sys.argv[1:])
dosages = args.dose
model = args.exp
samples = args.samples
outfile = args.out

#Run predict.py
if not os.path.isfile(outfile):
    predict_command = "python3 /home/wheelerlab3/MetaXcan/software/Predict.py "
    predict_command = predict_command + "--model_db_path "+model+" "
    predict_command = predict_command + "--model_db_snp_key varID "
    predict_command = predict_command + "--text_genotypes "+dosages+"nohead_hg38chr*.maf0.01.R20.8.dosage.txt "
    predict_command = predict_command + "--text_sample_ids "+samples+" "
    predict_command = predict_command + '--on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" '
    predict_command = predict_command + "--prediction_output "+outfile+".txt "
    predict_command = predict_command + "--prediction_summary_output "+outfile+"_summary.txt "
    #print(predict_command)
    os.system(predict_command)
