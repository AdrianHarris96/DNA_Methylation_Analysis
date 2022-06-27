#!/usr/bin/env python3

#Example input: ./automate_methylation.py -w /home/amharris/dna_methylation_analysis/

import subprocess as sp
import argparse as ap

parser = ap.ArgumentParser()
parser = ap.ArgumentParser(description='automate the methylation analysis')
parser.add_argument("-w", "--working", help="working directory with comparison sheets", required=True)
args = parser.parse_args()

comp_sheets = ['comp1.csv', 'comp2.csv', 'comp3.csv']
comp_list = ['K1_Low_K1_High', 'K2_Low_K2_High', 'K1_High_K2_High', 'K1_Low_K2_Low', 'K1_High_K2_Low', 'K1_Low_K2_High']


for sheet in comp_sheets:
	path = args.working + sheet
	for comp in comp_list:
		autorun = sp.call(['Rscript', 'methylation_analysis.R', path, '/local/projects-t3/XVMAS_Lab/Projects_2022/XVMAS_P09_methylation/01_analysis/kidneyTx_methylation/', '/home/amharris/dna_methylation_analysis/', '/local/projects-t3/XVMAS_Lab/Projects_2022/XVMAS_P09_methylation/01_analysis/kidneyTx_methylation/kidney_output/', comp])
