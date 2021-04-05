#Python script for running PredictAssociation.py
#Part of TWAS_for_protein Pipeline
#Created by Henry Wittich
'''This python script takes PrediXcan prediceted transcript levels and associates them with protein levels
for a given GTEx tissue, producing the following output:

tissue_type_proteinID_association.txt'''

import argparse
import sys
import os

cwd = os.getcwd() #current working directory

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Script to run PredictAssocation.py on every protein')
    parser.add_argument('--pred',
                        help='path to the predicted transcript levels',
                        required='True')
    parser.add_argument('--prot',
                        help='name of phenotype file (protein levels)',
                        required='True')
    parser.add_argument('--tissue',
                        help='name of the tissue model',
                        required='True')
    parser.add_argument('--out','-o',
                        help='path to output directory',
                        default=cwd)
    return parser.parse_args(args)

#retrieve command line arguments
args=check_arg(sys.argv[1:])
transcripts = args.pred
proteins = args.prot
protein_ID = proteins[:-4]
model = args.tissue
outfile = args.out

#Run PrediXcan
association_command = "python3 /home/wheelerlab3/MetaXcan/software/PrediXcanAssociation.py "
association_command = association_command + "--expression_file "+transcripts+" "
association_command = association_command + "--input_phenos_file "+outfile+"protein_levels/"+str(protein_ID)+".txt "
association_command = association_command + "--input_phenos_column "+str(protein_ID)+" "
association_command = association_command + "--output "+outfile+model+"/"+model+"_"+str(protein_ID)+"_association.txt "
#print("running association for "+str(protein_ID))
#print(association_command)
os.system(association_command)

