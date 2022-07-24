#! /usr/bin/env python3

import pandas as pd
import os
'''
df = pd.read_csv('/Users/adrianharris/Documents/dna_methylation_analysis/kidney_sample_sheet.csv') #Load in sample sheet (csv) with all samples

df = df.iloc[: , 1:] #Drop "unamed" sample 

new_sample_df = pd.DataFrame(columns = ['Basename', 'new_sample_id', 'array', 'eGFR1mo', 'eGFR12mo', 'eGFR24mo', 'donor_g']) #Creation of a new dataframe
for x in range(1, 3): #Number of sheets within excel file 
	sheet_df =  pd.read_excel('/Users/adrianharris/Desktop/kidney/Kidney Tx sample annotation-2.xlsx' , sheet_name=x, keep_default_na=False) #Updated sample sheet from Haseeb
	sheet_df.columns = [col.strip() for col in sheet_df.columns] 
	sample_id_column = sheet_df.columns[0][:-1] + str(x) #Important for pulling this column
	for index in sheet_df.index: #Creation of Basename and the new_sample_id
		basename =  str(sheet_df['Sentrix_ID'][index]) + '_' + str(sheet_df['Sentrix_Position'][index])
		new_sample_id = sheet_df[sample_id_column][index]
		array = sheet_df['Array_type'][index]
		gender = sheet_df['gender_donor'][index]
		if sheet_df['eGFR_1mo'][index] == '':
			print('Empty')
		elif float(sheet_df['eGFR_1mo'][index]) >= 45:
			e1 = 'High'
		else:
			e1 = 'Low'
		if sheet_df['eGFR_12mo'][index] == '':
			print('Empty')
		elif float(sheet_df['eGFR_12mo'][index]) >= 45:
			e12 = 'High'
		else:
			e12 = 'Low'
		if sheet_df['eGFR_24mo'][index] == '':
			print('Empty')
		elif float(sheet_df['eGFR_24mo'][index]) >= 45:
			e24 = 'High'
		else:
			e24 = 'Low'	
		if new_sample_id != '':
			new_sample_df.loc[len(new_sample_df)] = [basename, new_sample_id, array, e1, e12, e24, gender]

df = pd.merge(new_sample_df, df, on = 'Basename', how="right") #Right-outer join

#Drop old columns
df.drop(['array_type', 'eGFR_1month', 'eGFR_12month', 'eGFR_24month', 'donor_gender'], axis=1, inplace=True)

df.drop(columns='sample_id', inplace=True) #Drop the old sample_id column

#Rename new columns with old names
df = df.rename({'new_sample_id' : 'sample_id', 'array' : 'array_type', 'eGFR1mo': 'eGFR_1month', 'eGFR12mo' : 'eGFR_12month', 'eGFR24mo' : 'eGFR_24month', 'donor_g' : 'donor_gender'}, axis='columns') #rename the 'new_sample_id' to 'sample_id' for continuity

#Organize columns
df = df[['Basename', 'sample_id', 'sample_well', 'sample_plate', 'sentrix_ID', 'sentrix_position', 'time', 'array_type', 'donor_age', 'donor_gender', 'donor_race', 'recipient_age', 'recipient_gender', 'recipient_race', 'eGFR_1month', 'eGFR_12month', 'eGFR_24month', 'DGF_status', 'DCD_status']]

df.to_csv('/Users/adrianharris/Desktop/kidney_sample_sheet.csv', index = False)
'''
df = pd.read_csv('/Users/adrianharris/Documents/dna_methylation_analysis/kidney_sample_sheet.csv')
#df = df[df['eGFR_1month'].notna()] #Used to get paired samples with eGFR information
#print(df)
sample_dictionary = {}

for index in df.index:
	sample = df['sample_id'][index]
	if df['array_type'][index] == 'EPIC': #Only load EPIC samples 
		if sample[:-1] not in sample_dictionary.keys():
			sample_dictionary[sample[:-1]] = 1
		else:
			sample_dictionary[sample[:-1]] += 1 #Indicative of a paired sample if the final value is 2
	else:
		pass

sample_dictionary = dict(sorted(sample_dictionary.items()))

final_dictionary = {} #Load with keys that have a value-pair of 2 i.e. is paired
for key in sample_dictionary.keys():
    if sample_dictionary[key] == 2:
        final_dictionary[key] = 2

#print(final_dictionary)
#print(len(final_dictionary))

paired_df = pd.DataFrame(columns = df.columns) #Creation of paired_df using the final_dictionary and the df

for index in df.index:
	sample = df['sample_id'][index]
	if sample[:-1] in final_dictionary.keys():
		paired_df.loc[len(paired_df)] = list(df.loc[index,])

paired_df.to_csv('/Users/adrianharris/Desktop/paired_kidney_sample_sheet.csv', index = False)