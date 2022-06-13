#! /usr/bin/env python3

#Example input: ./sort_array_type.py -x /Users/adrianharris/Downloads/Kidney\ Tx\ sample\ annotation.xlsx

import pandas as pd
import argparse as ap
import os
import shutil 

parser = ap.ArgumentParser()
parser = ap.ArgumentParser(description='Calculate effect dose for each UKBB individual using the matrix')
parser.add_argument("-x", "--excel", help="excel sample file", required=True)
parser.add_argument("-w", "--working", help="working directory", default='/Users/adrianharris/Desktop/kidney/')
parser.add_argument("-o", "--output", help="output csv filename", default='methylation.csv')
args = parser.parse_args()

arrayDict = {} #Dictionary intended for moving files later in script

df = pd.DataFrame(columns=['sample_id', 'sample_well', 'sample_plate', 'sentrix_ID', 'sentrix_position', 'time', 'array_type', 'gender_donor', 'gender_recipient'])
for x in range(1, 5): #Number of sheets within excel file 
	sheet_df =  pd.read_excel(args.excel, sheet_name=x, keep_default_na=True)
	sheet_df.columns = [col.strip() for col in sheet_df.columns] #Stripping whitespace in column names
	if x == 4: #Used for the extraction of the sample_id from the excel file 
		x = 12
	sample_id = 'SampleID_K' + str(x)
	for index in sheet_df.index: #iterate through columns
		row = sheet_df.loc[index, [sample_id, 'Sample_Well', 'Sample_Plate', 'Sentrix_ID', 'Sentrix_Position', 'Time', 'Array_type', 'gender_donor', 'gender_recipient']]
		row = list(row)
		sentrix = str(sheet_df['Sentrix_ID'][index]) + '_' + str(sheet_df['Sentrix_Position'][index])
		array = sheet_df['Array_type'][index]
		if sentrix in arrayDict.keys():
			pass
			#print('duplicate')
		else:
			arrayDict[sentrix] = array
		df.loc[len(df) + 1,] = row #adding to df

#print(arrayDict)

df = df[df['time'].notna()] #drop if time is empty
#print(df)

#Appending to either 450K or EPIC list 
list450K = []
listEPIC = []
for key in arrayDict.keys():
	if arrayDict[key] == '450K':
		list450K.append(key)
	elif arrayDict[key] == 'EPIC':
		listEPIC.append(key)
	else:
		pass

#print(len(list450K))
#print(len(listEPIC))

#Creation of new sample sheet (CSV) using output argument
os.chdir(args.working)
df.to_csv(args.output)

#Determine all files
f = []
for filename in os.listdir(args.working):
    f.append(filename)

os.makedirs('450k_array')
os.makedirs('EPIC_array')

#Moving to corresponding directory based on array_type
for item in list450K:
	for file in f:
		if item in file:
			src = args.working + file
			dst = args.working + '450k_array/' + file
			shutil.copyfile(src, dst)

for item in listEPIC:
	for file in f:
		if item in file:
			src = args.working + file
			dst = args.working + 'EPIC_array/' + file
			shutil.copyfile(src, dst)

