#! /usr/bin/env python3

#Example input: ./sort_array_type.py -x /Users/adrianharris/Desktop/kidney/kidney\ Tx\ sample\ annotation.xlsx -s 1 -e 5

import pandas as pd
import argparse as ap
import os
import shutil 

parser = ap.ArgumentParser()
parser = ap.ArgumentParser(description='Sort raw files into respective folder according to array type')
parser.add_argument("-x", "--excel", help="excel sample file", required=True)
parser.add_argument("-w", "--working", help="working directory", default='/Users/adrianharris/Desktop/kidney/')
parser.add_argument("-s", "--start", help="starting sheet(#) in excel file", type=int, required=True)
parser.add_argument("-e", "--end", help="ending sheet(#) in excel file", type=int, required=True)
args = parser.parse_args()

arrayDict = {} #Dictionary intended for moving files later in script

for x in range(args.start, args.end): #Iterate through important sheets in the excel file 
	sheet_df =  pd.read_excel(args.excel, sheet_name=x, keep_default_na=True)
	sheet_df.columns = [col.strip() for col in sheet_df.columns] #Stripping whitespace in column names
	for index in sheet_df.index:
		sentrix = str(sheet_df['Sentrix_ID'][index]) + '_' + str(sheet_df['Sentrix_Position'][index])
		array = sheet_df['Array_type'][index]
		if sentrix in arrayDict.keys():
			pass #avoid loading duplicates with this line
		else:
			arrayDict[sentrix] = array

#print(arrayDict)

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

#Move to working directory and sort files into respective folders 
os.chdir(args.working)

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

