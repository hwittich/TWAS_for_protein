#Code for running PrediXcan's association script using observed expression data instead of
#predicted
#Adapted from by TWAS for protein pipeline

import argparse
import sys
import subprocess
import os

cwd = os.getcwd() #current working directory

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Script to run PrediXcan for every protein in every tissue.')
    parser.add_argument('--prot',
                        help='path to the protein levels',
                        required='True')
    parser.add_argument('--tissue',
                        help='tissue',
                        required='True')                        
    parser.add_argument('--samples',
                        help='path to the list of samples',
                        required='True')
    parser.add_argument('--out','-o',
                        help='path to output directory',
                        default=cwd)
    return parser.parse_args(args)

#retrieve command line arguments
args=check_arg(sys.argv[1:])
proteins = cwd+"/"+args.prot
tissue = args.tissue
samples = cwd+"/"+args.samples
outfile = cwd+"/"+args.out
if not (outfile[-1] == "/"):
    outfile = outfile + "/"
if not os.path.isdir(outfile):
    os.system("mkdir -p "+outfile)
    
###Reformatting some of the input files
#Creating samples.txt file
protein_matrix = open(proteins,'r')
lines = protein_matrix.readlines() #get list of line in the file
sample_IDs = lines[0].split(" ") #header of each column is a sample ID #space separated
if not os.path.isfile(samples):
    samples_file = open(samples,'w')
    for ID in sample_IDs[1:]:
        #print(ID)
        samples_file.write(str(ID.rstrip())+"\t"+str(ID.rstrip())+"\n")
    samples_file.close()
    
#Divide up protein matrix by protein
#Set up output directory
if not os.path.isdir(outfile+"protein_levels/"):
    os.system("mkdir "+outfile+"protein_levels")
SL_IDs={} #Store the protein IDs for use later in pipeline
multigene_aptamers=set() #Store the multigene aptamers for use later in pipeline
for line in lines[1:]: #each row pertains to 1 protein. #skip first line with headers
    #Make a temporary file where column 1 is sample_IDs and column two is protein levels for current protein
    protein_levels = line.split(" ") #split line around separator
    joint_ID = protein_levels[0].split("_") #Joint ID format <somascan (SL) ID>_<ensembl ID> #Take the SL ID
    SL_ID = joint_ID[0].rstrip()
    if SL_ID in SL_IDs.keys():
        multigene_aptamers.add(str(SL_ID)+"_"+str(SL_IDs[SL_ID])) #Store the original joint ID in the multigene aptamers set
        multigene_aptamers.add(str(SL_ID)+"_"+str(joint_ID[1].rstrip()))
    else: #Otherwise, new aptamer ID, store ENSG ID as value
        SL_IDs[SL_ID]=joint_ID[1].rstrip()
    protein_file = open(outfile+"protein_levels/"+SL_ID+"_"+tissue+".txt",'w')
    #Line 1, write headers
    protein_file.write("FID\tIID\t"+str(SL_ID)+"\n")
    #Move through list of sample IDs (same length as list of protein levels) and write the sample ID followed by the protein level
    for i in range(1,len(sample_IDs)):
        protein_file.write(str(sample_IDs[i]).rstrip()+"\t"+str(sample_IDs[i]).rstrip()+"\t"+str(protein_levels[i]).rstrip()+"\n")
    #close protein file
    protein_file.close()

#Write list of multigene aptamers to output file
multigene_aptamers_file=open(outfile+"multigene_aptamers.txt",'w')
multigene_aptamers_file.write("joint_ID\n")
for joint_ID in multigene_aptamers:
    multigene_aptamers_file.write(str(joint_ID)+"\n")
multigene_aptamers_file.close()

###Running PrediXcan
#Set up output directory
if not os.path.isdir(outfile+tissue):
    os.system("mkdir "+outfile+tissue)

##Perform PredictAssociation.py
#Split up job over 15 cores
cores=15
processes=set()
for SL_ID in SL_IDs: #Loop through every protein
    association_command = subprocess.Popen(["python3","code/protein_association.py","--pred",outfile+tissue+"_obs_expression.txt","--phen",SL_ID+"_"+tissue+".txt","--prot",SL_ID,"--tissue",tissue,"--out",outfile])
    processes.add(association_command)
    if len(processes) >= cores: #Make sure number of processes isn't greater than number of cores
        os.wait() #Wait for a process to finish before starting a new one
        processes.difference_update([p for p in processes if p.poll() is not None])

#Make sure all association commands have finished before continuing
for p in processes:
    if p.poll() is None:
        p.wait()        
   
###Analyze the results of PrediXcan
##Concatenate PrediXcan output files
#Loop through the association files, then loop through each line and write to output
#Don't include headers
every_protein_association = open(outfile+tissue+"/obs_expr_"+tissue+"_every_protein_association.txt",'w')
every_protein_association.write("protein\tgene\teffect\tse\tzscore\tpvalue\tn_samples\tstatus\n")
for protein in SL_IDs:
    protein_association = open(outfile+tissue+"/"+tissue+"_"+str(protein)+"_association.txt",'r')
    results = protein_association.readlines()
    for line in results[1:]:
        every_protein_association.write(str(protein)+"\t"+str(line))
    protein_association.close()
every_protein_association.close()
protein_matrix.close()
