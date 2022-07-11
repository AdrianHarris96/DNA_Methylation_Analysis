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
	sent = int(sheet_df['Sentrix_ID'][index])
	basename = str(sent) + '_' + str(sheet_df['Sentrix_Position'][index])
	row = sheet_df.loc[index, ['Sample_Name', 'Donor type', 'Donor age', 'Donor gender', 'Donor race', 'Sample_Well', 'Sample_Plate', 'Sample_Group', 'Sentrix_ID', 'array type', 'Sentrix_Position']]
	row = list(row)
	row.insert(0, basename)
	df.loc[len(df) + 1,] = row #adding to df

df.drop_duplicates(inplace=True)

#Create new column for L1 and L2 
df['collection'] = 'NA'
for index in df.index:
	if 'L1' in df['sample_name'][index]:
		df['collection'][index] = 'L1'
	else:
		df['collection'][index] = 'L2'

#Loading donor_type_dict
donor_type_dict = {}
for index in df.index:
	if df['donor_type'][index] == 'DD' or df['donor_type'][index] == 'DCD' or df['donor_type'][index] == 'LD':
		sample = df['sample_name'][index]
		if ' ' in sample:
			key = sample.split(' ')[0]
			key = key.strip()
		if 'L' in sample:
			key = sample.split('L')[0]
			key = key.strip()
		donor_type_dict[key] = [df['donor_type'][index], df['donor_age'][index], df['donor_gender'][index], df['donor_race'][index]]
	else:
		pass

#Filling attributes in the donor_type column
for index in df.index:
	if df['collection'][index] == 'L2':
		key = df['sample_name'][index][0:4]
		if key in donor_type_dict.keys():
			df['donor_type'][index] = donor_type_dict[key][0]
			df['donor_age'][index] = donor_type_dict[key][1]
			df['donor_gender'][index] = donor_type_dict[key][2]
			df['donor_race'][index] = donor_type_dict[key][3]
		else:
			pass
	else:
		pass
os.chdir(args.out_dir) #Move to output_directory

#print(df)
df.to_csv(args.output)
