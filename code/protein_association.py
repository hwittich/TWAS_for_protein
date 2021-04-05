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
                        help='path to the protein levels',
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
model = args.tissue
outfile = args.out

#Loop through the lines of the protein file
protein_matrix = open(proteins,'r')
lines = protein_matrix.readlines() #get list of line in the file
sample_IDs = lines[0].split(" ") #header of each column is a sample ID #space separated

#Set up output directory
if not os.path.isdir(outfile+"protein_levels/"):
    os.system("mkdir "+outfile+"protein_levels")
for row in lines[1:]: #each row pertains to 1 protein. #skip first line with headers
    #Make a temporary file where column 1 is sample IDs and column two is protein levels for current protein
    protein_levels = row.split(" ") #split line around separator
    protein_ID = str(protein_levels[0].split("_")[0]).rstrip() #Joint ID format <somascan (SL) ID>_<ensembl ID> #Take the SL ID
    protein_file = open(outfile+"/protein_levels/"+protein_ID+".txt",'w') #open file to write transpose protein levels
    #Line 1, write headers
    protein_file.write("FID\tIID\t"+str(protein_ID)+"\n")
    #Move through list of sample IDs (same length as list of protein levels) and write the sample ID followed by the protein level
    for i in range(1,len(sample_IDs)):
        protein_file.write(str(sample_IDs[i]).rstrip()+"\t"+str(sample_IDs[i]).rstrip()+"\t"+str(protein_levels[i]).rstrip()+"\n")
    #close protein file
    protein_file.close()

    #Now use this protein expression file to run prediXcan!
    association_command = "python3 /home/wheelerlab3/MetaXcan/software/PrediXcanAssociation.py "
    association_command = association_command + "--expression_file "+pred+" "
    association_command = association_command + "--input_phenos_file "+outfile+"protein_levels/"+str(protein_ID)+".txt "
    association_command = association_command + "--input_phenos_column "+str(protein_ID)+" "
    association_command = association_command + "--output "+outfile+model+"_"+str(protein_ID)+"_association.txt "
    #print("running association for "+str(protein_ID))
    #print(association_command)
    os.system(association_command)

