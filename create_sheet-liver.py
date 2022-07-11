#! /usr/bin/env python3

#Example input: ./create_sheet-liver.py -x /Users/adrianharris/Desktop/liver/liver_samples_annotation.xlsx -o liver_sample_sheet.csv -d /Users/adrianharris/Documents/dna_methylation_analysis/

import pandas as pd
import argparse as ap
import os

parser = ap.ArgumentParser()
parser = ap.ArgumentParser(description='Create sample sheet')
parser.add_argument("-x", "--excel", help="excel sample file", required=True)
parser.add_argument("-d", "--out_dir", help="output directory", default='/Users/adrianharris/Desktop/kidney/')
parser.add_argument("-o", "--output", help="output csv filename", default='methylation.csv')
args = parser.parse_args()

df = pd.DataFrame(columns=['Basename', 'sample_name', 'donor_type', 'donor_age', 'donor_gender', 'donor_race', 'sample_well', 'sample_plate', 'sample_group', 'sentrix_ID', 'array_type', 'sentrix_position'])

sheet_df =  pd.read_excel(args.excel, sheet_name=0, keep_default_na=True)
sheet_df.columns = [col.strip() for col in sheet_df.columns] #Stripping whitespace in column names
for index in sheet_df.index: #iterate through columns
	basename = str(sheet_df['Sentrix_ID'][index])[:-2] + '_' + str(sheet_df['Sentrix_Position'][index])
	row = sheet_df.loc[index, ['Sample_Name', 'Donor type', 'Donor age', 'Donor gender', 'Donor race', 'Sample_Well', 'Sample_Plate', 'Sample_Group', 'Sentrix_ID', 'array type', 'Sentrix_Position']]
	row = list(row)
	row.insert(0, basename)
	df.loc[len(df) + 1,] = row #adding to df

df.drop_duplicates(inplace=True)

os.chdir(args.out_dir) #Move to output_directory

df.to_csv(args.output)
