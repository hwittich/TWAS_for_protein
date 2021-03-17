#Python Wrapper for TWAS_for_protein Pipeline
#Created by Henry Wittich
'''This python script takes genotype dosage data and protein expression data
from a population in the TOPMed cohort, along with a gene expression model
and runs PrediXcan to  conduct association tests between the genotypes of the
individuals and the expression levels of every protein, producing the following
output file:

population_tissue_type_protein_association.txt

For usage, type from the command line:
python protein_association.py'''

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
    parser = argparse.ArgumentParser(description='Script to run PrediXcan on all protein')
    parser.add_argument('--dosages',
                        help='path to the dosage file',
                        required='True')
    parser.add_argument('--exp',
                        help='path to the expression models.',
                        required='True')
    parser.add_argument('--proteins',
                        help='path to the protein levels',
                        required='True')
    parser.add_argument('--samples',
                        help='path to the list of samples',
                        required='True')
    parser.add_argument('--out','-o',
                        help='file to output as',
                        default=cwd)
    parser.add_argument('--software',
                        help='path to the directory containing the PrediXcan scripts',
                        required='True')
    return parser.parse_args(args)

#retrieve command line arguments
args=check_arg(sys.argv[1:])
dosages = cwd+"/"+args.dosages
models = cwd+"/"+args.exp
proteins = cwd+"/"+args.proteins
samples = cwd+"/"+args.samples
outfile = cwd+"/"+args.out
scripts = args.software

#Loop through the lines of the protein file
protein_matrix = open(proteins,'r')
lines = protein_matrix.readlines() #get list of line in the file
sample_IDs = lines[0].split(" ") #header of each column is a sample ID #space separated
#Making samples.txt file
samples_file = open(samples,'w')
for ID in sample_IDs[1:]:
    print(ID)
    samples_file.write(str(ID.rstrip())+"\t"+str(ID.rstrip())+"\n")
samples_file.close()

#Loop through the dosage files and remove the headers
for chromosome in range(1,23):
    dosage_in = open(dosages+"/hg38chr"+str(chromosome)+".maf0.01.R20.8.dosage.txt","r")
    dosage_out = open(dosages+"/nohead_hg38chr"+str(chromosome)+".maf0.01.R20.8.dosage.txt","w")
    for line in dosage_in.readlines()[1:]: #Write every line but first to a new output file
        dosage_out.write(line)
    dosage_in.close()
    dosage_out.close()

#First, run the predict.py code
predict_command = "python3 "+scripts+"/Predict.py "
predict_command = predict_command + "--model_db_path "+models+" "
predict_command = predict_command + "--model_db_snp_key varID "
predict_command = predict_command + "--text_genotypes "+dosages+"/nohead_hg38chr*.maf0.01.R20.8.dosage.txt "
predict_command = predict_command + "--text_sample_ids "+samples+" "
predict_command = predict_command + '--on_the_fly_mapping METADATA "chr{}_{}_{}_{}_b38" '
predict_command = predict_command + "--prediction_output "+outfile+"/Whole_Blood_predict.txt "
predict_command = predict_command + "--prediction_summary_output "+outfile+"/Whole_Blood_predict_summary.txt "
#print(predict_command)
os.system(predict_command)


for row in lines[1:]: #each row pertains to 1 protein. #skip first line with headers
    #Make a temporary file where column 1 is sample IDs and column two is protein levels for current protein
    protein_levels = row.split(" ") #split line around separator
    protein_ID = str(protein_levels[0].split("_")[1]).rstrip() #Joint ID format <soma ID>_<ensembl ID> #Take the ensembl ID
    protein_file = open(outfile+"/protein_levels/"+protein_ID+".txt",'w') #open file to write transpose protein levels
    #Line 1, write headers
    protein_file.write("FID\tIID\t"+str(protein_ID)+"\n")
    #Move through list of sample IDs (same length as list of protein levels) and write the sample ID followed by the protein level
    for i in range(1,len(sample_IDs)):
        protein_file.write(str(sample_IDs[i]).rstrip()+"\t"+str(sample_IDs[i]).rstrip()+"\t"+str(protein_levels[i]).rstrip()+"\n")

    #Now use this protein expression file to run prediXcan!
    association_command = "python3 "+scripts+"/PrediXcanAssociation.py "
    association_command = association_command + "--expression_file "+outfile+"/Whole_Blood_predict.txt "
    association_command = association_command + "--input_phenos_file "+outfile+"/protein_levels/"+protein_ID+".txt "
    association_command = association_command + "--input_phenos_column "+str(protein_ID)+" "
    association_command = association_command + "--output "+outfile+"/Whole_Blood_"+str(protein_ID)+"_association.txt "
    os.system(association_command)

    #close protein file
    protein_file.close()
