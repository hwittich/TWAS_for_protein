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

import argparse
import sys
import os

cwd = os.getcwd() #current working directory

def check_arg(args=None):
    #Thank you Ryan Schubert for some of this code
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
    return parser.parse_args(args)

#retrieve command line arguments
args=check_arg(sys.argv[1:])
dosages = args.dosages
models = args.expression_models
proteins = args.proteins
samples = args.samples
outfile = args.out
