#Python Wrapper for TWAS_for_protein Pipeline
#Created by Henry Wittich
'''This python script takes genotype dosage data and protein expression data
from a population in the TOPMed cohort, along with a gene expression model
and runs PrediXcan to  conduct association tests between the genotypes of the
individuals and the expression levels of every protein. 

For usage, type from the command line:
python3 protein_TWAS.py'''

#python3 protein_association.py --dosages home/hwittich/mount/TWAS_for_protein/TOPMed_ALL_data/dosages --exp /home/hwittich/mount/TWAS_for_protein/TOPMed_Test_Data/en_Whole_Blood.db --proteins /home/hwittich/mount/TWAS_for_protein/TOPMed_ALL_data/Proteome_TOPMed_ALL_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids_sorted.txt --samples /home/hwittich/mount/TWAS_for_protein/TOPMed_ALL_data/samples.txt --out /home/hwittich/mount/TWAS_for_protein/output/3-5-2021
#python3 code/protein_association.py --dosages TOPMed_ALL_data/dosages --exp TOPMed_Test_Data/en_Whole_Blood.db --proteins TOPMed_ALL_data/Proteome_TOPMed_ALL_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids_sorted.txt --samples TOPMed_ALL_data/samples.txt --out output/3-5-2021 --software /home/wheelerlab3/MetaXcan/software
#nohup python3 code/protein_association.py --dosages TOPMed_ALL_data/dosages --exp TOPMed_Test_Data/en_Whole_Blood.db --proteins TOPMed_ALL_data/Proteome_TOPMed_ALL_ln_adjAgeSex_mean_rank-inverse_adj10PCair_MEQTL_reformatted_ids_sorted.txt --samples TOPMed_ALL_data/samples.txt --out output/3-12-2021 --software /home/wheelerlab3/MetaXcan/software &
#PredictAssociation Command
#python3 /home/wheelerlab3/MetaXcan/software/PrediXcanAssociation.py --expression_file /home/henry/TWAS_for_protein/output/2-3-2021/Whole_Blood_predict.txt --input_phenos_file /home/henry/TWAS_for_protein/output/3-12-2021/protein_levels/ENSG00000111674.9.txt --input_phenos_column ENSG00000111674.9 --output /home/henry/TWAS_for_protein/output/3-12-2021/Whole_Blood_association.txt

import argparse
import sys
import os

cwd = os.getcwd() #current working directory

def check_arg(args=None):
    parser = argparse.ArgumentParser(description='Script to run PrediXcan for every protein in every tissue.')
    parser.add_argument('--dose',
                        help='path to the dosage files',
                        required='True')
    parser.add_argument('--exp',
                        help='path to the expression models',
                        required='True')
    parser.add_argument('--prot',
                        help='path to the protein levels',
                        required='True')
    parser.add_argument('--samples',
                        help='path to the list of samples',
                        required='True')
    parser.add_argument('--genes',
                        help='path to gene annotation file',
                        required='True')
    parser.add_argument('--aptamers',
                        help='path to protein annotation file',
                        required='True')
    parser.add_argument('--out','-o',
                        help='path to output directory',
                        default=cwd)
    return parser.parse_args(args)

#retrieve command line arguments
args=check_arg(sys.argv[1:])
dosages = cwd+"/"+args.dose
if not (dosages[-1] == "/"):
    dosages = dosages + "/"
models = cwd+"/"+args.exp
if not (models[-1] == "/"):
    models = models + "/"
proteins = cwd+"/"+args.prot
samples = cwd+"/"+args.samples
gene_anno = cwd+"/"+args.genes
protein_anno = cwd+"/"+args.aptamers
outfile = cwd+"/"+args.out
if not (outfile[-1] == "/"):
    outfile = outfile + "/"
if not os.path.isdir(outfile):
    os.system("mkdir -p "+outfile)

#print(dosages)
#print(models)
#print(proteins)
#print(samples)
#print(outfile)

###Reformatting some of the input files
#Creating samples.txt file
if not os.path.isfile(samples):
    protein_matrix = open(proteins,'r')
    lines = protein_matrix.readlines() #get list of line in the file
    protein_matrix.close()
    sample_IDs = lines[0].split(" ") #header of each column is a sample ID #space separated


    samples_file = open(samples,'w')
    for ID in sample_IDs[1:]:
        #print(ID)
        samples_file.write(str(ID.rstrip())+"\t"+str(ID.rstrip())+"\n")
    samples_file.close()

#Remove the headers from the dosage files
for chromosome in range(1,23):
    dosage_in = open(dosages+"/hg38chr"+str(chromosome)+".maf0.01.R20.8.dosage.txt","r")
    dosage_out = open(dosages+"/nohead_hg38chr"+str(chromosome)+".maf0.01.R20.8.dosage.txt","w")
    for line in dosage_in.readlines()[1:]: #Write every line but first to a new output file
        dosage_out.write(line)
    dosage_in.close()
    dosage_out.close()

###Running PrediXcan
#Loop through every tissue model
for filename in os.listdir(models):
    if filename[-3:] == ".db":
        ##Perform predict.py
        os.system("python3 code/transcript_prediction.py --dosages "+dosages+" --exp "+models+filename+" --samples "+samples+" --out "+outfile+filename[:-3]+"_predict")

        ##Perform PredictAssociation.py
        #Split up job over 20 cores
        #Loop through the lines of the protein matrix and split up file into 20 partitions
        #First set up output directory
        partition_dir = outfile+"protein_matrix_partitions/"
        if not os.path.isdir(partition_dir):
            os.system("mkdir "+partition_dir)
        #Create first partition
        partition=1 #current partition of the proteins
        partition_path = partition_dir+"proteome_"+str(partition)+".txt"
        protein_partition = open(partition_path,'w')
        protein_partition.write(lines[0]) #Write header line of protein matrix to file

        count=0 #number of proteins in each file
        #Reset count once it reaches (total # proteins)/(19 cores) #remainder goes on final core
        SL_IDs={} #Store the protein IDs for use later in pipeline
        multigene_aptamers=set() #Store the multigene aptamers for use later in pipeline
        for line in lines[1:]: #first line is the sample IDs
            count += 1 #increment count
            joint_ID = line.split(" ")[0].split("_")
            SL_ID = joint_ID[0].rstrip()
            if SL_ID in SL_ID.keys():
                multigene_aptamers.add(str(SL_ID)+"_"+str(SL_IDs[SL_ID])) #Store the original joint ID in the multigene aptamers set
                multigene_aptamers.add(str(SL_ID)+"_"+str(joint_ID[1].rstrip()))
            else: #Otherwise, new aptamer ID, store ENSG ID as value
                SL_IDs[SL_ID]=joint_ID[1].rstrip()
                
            if count == (len(lines)-1)//20:
                #Close output file
                protein_partition.close()

                #Running PredictAssociation.py on completed partition
                os.system("nohup python3 code/protein_association.py --pred "+outfile+filename[:-3]+"_predict.txt --prot "+partition_path+" --out "+outfile+" --tissue "+filename[:-3]+" &")

                #Start next partition
                partition+=1
                partition_path = partition_dir+"proteome_"+str(partition)+".txt"
                protein_partition = open(partition_path,'w')
                protein_partition.write(lines[0])
                count=1 #reset count
            protein_partition.write(line)
        protein_partition.close()
        #Run predictAssociation.py on final partition
        os.system("nohup python3 code/protein_association.py --pred "+outfile+filename[:-3]+"_predict.txt --prot "+partition_path+" --out "+outfile+" --tissue "+filename[:-3]+" &")

        #Write list of multigene aptamers to output file
        multigene_aptamers_file=open(outfile+"protein_matrix_partitions/multigene_aptamers.txt",'w')
        multigene_aptamers_file.write("joint_ID\n")
        for joint_ID in multigene_aptamers:
            multigene_aptamers_file.write(str(joint_ID)+"\n")
        multigene_aptamers_file.close()
        
        ###Analyze the results of PrediXcan
        ##Concatenate PrediXcan output files
        #Loop through the association files, then loop through each line and write to output
        #Don't include headers
        every_protein_association = open(outfile+filename[:-3]+"_every_protein_association.txt",'w')
        every_protein_association.write("protein\tgene\teffect\tse\tzscore\tpvalue\tn_samples\tstatus\n")
        for protein in SL_IDs:
            protein_association = open(outfile+filename[:-3]+"_"+str(protein)+"_association.txt",'r')
            results = protein_association.readlines()
            for line in results[1:]:
                every_protein_association.write(str(protein)+"\t"+str(line))
            protein_association.close()
        every_protein_association.close()

        ##Plot in R
        #Set up output directory
        if not os.path.isdir(outfile+"significant_pairs_plots/"):
            os.system("mkdir "+outfile+"significant_pairs_plots/")
        os.system("Rscript code/plot_significant_pairs.R "+outfile+filename[:-3]+"_every_protein_association.txt "+filename[:-3]+" "+gene_anno+" "+protein_anno+" "+outfile+"protein_matrix_partitions/multigene_aptamers.txt "+outfile+"significant_pairs_plots/")
        
