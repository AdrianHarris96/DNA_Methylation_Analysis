#! /usr/bin/env python3

#Example input: ./create_sheet-kidney.py -x /Users/adrianharris/Desktop/kidney/Kidney\ Tx\ sample\ annotation-3.xlsx -s 1 -e 3 -o /Users/adrianharris/Documents/dna_methylation_analysis/

import pandas as pd
import argparse as ap
import os

parser = ap.ArgumentParser()
parser = ap.ArgumentParser(description='Create sample sheet and paired for kidney')
parser.add_argument("-x", "--excel", help="excel sample file", required=True)
parser.add_argument("-o", "--out_dir", help="output directory", default='/Users/adrianharris/Desktop/kidney/')
parser.add_argument("-s", "--start", help="starting sheet(#) in excel file", type=int, required=True)
parser.add_argument("-e", "--end", help="ending sheet(#) in excel file", type=int, required=True)
parser.add_argument("-t1", "--time1", help="time point 1 to extract", default='K1')
parser.add_argument("-t2", "--time2", help="time point 2 to extract", default='K2')
args = parser.parse_args()

df = pd.DataFrame(columns=['Basename', 'sample_id', 'sample_well', 'sample_plate', 'sentrix_ID', 'sentrix_position', 'time', 'array_type', 'donor_gender', 'recipient_gender', 'eGFR_1month', 'eGFR_12month', 'eGFR_24month', 'ICT', 'DGF_status', 'DCD_status'])
for x in range(args.start, args.end): #Number of sheets within excel file 
	sheet_df =  pd.read_excel(args.excel , sheet_name=x, keep_default_na=False) #Updated sample sheet from Haseeb
	sheet_df.columns = [col.strip() for col in sheet_df.columns] 
	sample_id_column = sheet_df.columns[0][:-1] + str(x) #Important for pulling this column
	for index in sheet_df.index: 
		basename =  str(sheet_df['Sentrix_ID'][index]) + '_' + str(sheet_df['Sentrix_Position'][index]) #Creation of Basename 
		if (sheet_df['Time'][index] == args.time1) or (sheet_df['Time'][index] == args.time2): #Only append to df if either time point is present 
			basename = str(sheet_df['Sentrix_ID'][index]) + '_' + str(sheet_df['Sentrix_Position'][index])
			try:
				row = sheet_df.loc[index, [sample_id_column, 'Sample_Well', 'Sample_Plate', 'Sentrix_ID', 'Sentrix_Position', 'Time', 'Array_type', 'gender_donor', 'gender_recipient', 'eGFR_1mo', 'eGFR_12mo', 'eGFR_24mo', 'Ischemia Time (min)', 'DGF=1, non-DGF=0', 'DCD(yes/no)']]
				row = list(row)
				row.insert(0, basename)
			except:
				row = sheet_df.loc[index, [sample_id_column, 'Sample_Well', 'Sample_Plate', 'Sentrix_ID', 'Sentrix_Position', 'Time', 'Array_type', 'gender_donor', 'gender_recipient', 'eGFR_1mo', 'eGFR_12mo', 'eGFR_24mo']]
				row = list(row)
				row.extend(['NaN', 'NaN', 'NaN'])
				row.insert(0, basename)
			finally:
				df.loc[len(df) + 1,] = row #adding to df
		else:
			pass

del sheet_df

#Changing the value to a categorical, high-low
monthList = ['eGFR_1month', 'eGFR_12month', 'eGFR_24month']
for month in monthList:
	for index in df.index:
		try:
			if df[month][index] < 45:
				df[month][index] = 'Low'
			elif df[month][index] >= 45:
				df[month][index] = 'High'
		except:
			df.loc[index, month] = ''

df.to_csv('kidney_sample_sheet.csv', index = False)

sample_dictionary = {}
for index in df.index:
	sample = df['sample_id'][index]
	if sample[:-1] not in sample_dictionary.keys():
		sample_dictionary[sample[:-1]] = 1
	else:
		sample_dictionary[sample[:-1]] += 1 #Indicative of a paired sample if the final value is 2
	

sample_dictionary = dict(sorted(sample_dictionary.items()))

final_dictionary = {} #Load with keys that have a value-pair of 2 i.e. is paired
for key in sample_dictionary.keys():
    if sample_dictionary[key] == 2:
        final_dictionary[key] = 2

paired_df = pd.DataFrame(columns = df.columns) #Creation of paired_df using the final_dictionary and the df

#Filling the paired dataframe 
for index in df.index:
	sample = df['sample_id'][index]
	if sample[:-1] in final_dictionary.keys():
		paired_df.loc[len(paired_df)] = list(df.loc[index,])

del df

old_df = pd.DataFrame(columns=['Basename', 'donor_age', 'donor_race', 'recipient_age', 'recipient_race'])
#Merge with older dataframe
for x in range(0, 2):
	sheet_df =  pd.read_excel('/Users/adrianharris/Desktop/kidney/samples annotation.xlsx' , sheet_name=x, keep_default_na=False) #Updated sample sheet from Haseeb
	sheet_df.columns = [col.strip() for col in sheet_df.columns] 
	for index in sheet_df.index: 
		basename =  str(sheet_df['Sentrix_ID'][index]) + '_' + str(sheet_df['Sentrix_Position'][index])
		row = sheet_df.loc[index, ['Donor Age', 'Donor Race', 'Recipient Age', 'Recipient Race']]
		row = list(row)
		row.insert(0, basename)
		old_df.loc[len(old_df) + 1,] = row #adding to df

#Merge the paired dataframe with older dataframe 
paired_df = pd.merge(paired_df, old_df, on = 'Basename') 

#Ensure the K1 Covariates match that of K2
covar_dictionary = {}
for index in paired_df.index:
	sample = paired_df['sample_id'][index]
	sample = sample[:-1]
	if paired_df['time'][index] == 'K1':
		row = paired_df.loc[index, ['donor_gender', 'recipient_gender', 'eGFR_1month', 'eGFR_12month', 'eGFR_24month', 'ICT', 'DGF_status', 'DCD_status', 'donor_age', 'donor_race', 'recipient_age', 'recipient_race']]
		row = list(row)
		covar_dictionary[sample] = row
	else:
		paired_df.loc[index, ['donor_gender', 'recipient_gender', 'eGFR_1month', 'eGFR_12month', 'eGFR_24month', 'ICT', 'DGF_status', 'DCD_status', 'donor_age', 'donor_race', 'recipient_age', 'recipient_race']] = covar_dictionary[sample]

paired_df.loc[paired_df['sample_id'] == '67_K2', 'array_type'] = 'EPIC' #This sample is labeled incorrectly

df_450k = paired_df.loc[paired_df['array_type'] == '450K']
df_EPIC = paired_df.loc[paired_df['array_type'] == 'EPIC']
os.chdir(args.out_dir) #Move to output_directory

paired_df.to_csv('paired_kidney_sample_sheet.csv', index = False)
df_450k.to_csv('paired_450kkidney_sample_sheet.csv', index = False)
df_EPIC.to_csv('paired_EPICkidney_sample_sheet.csv', index = False)

quit()
