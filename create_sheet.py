#! /usr/bin/env python3

#Example input: ./create_sheet.py -x /Users/adrianharris/Desktop/kidney/Kidney\ Tx\ sample\ annotation.xlsx -s 1 -e 5 -o comp1.csv -m 

import pandas as pd
import argparse as ap
import os

parser = ap.ArgumentParser()
parser = ap.ArgumentParser(description='Create sample sheet')
parser.add_argument("-x", "--excel", help="excel sample file", required=True)
parser.add_argument("-w", "--working", help="working directory", default='/Users/adrianharris/Desktop/kidney/')
parser.add_argument("-o", "--output", help="output csv filename", default='methylation.csv')
parser.add_argument("-s", "--start", help="starting sheet(#) in excel file", type=int, required=True)
parser.add_argument("-e", "--end", help="ending sheet(#) in excel file", type=int, required=True)
parser.add_argument("-t1", "--time1", help="time point 1 to extract", default='K1')
parser.add_argument("-t2", "--time2", help="time point 2 to extract", default='K2')
parser.add_argument("-m", "--month", help="eGFR month", required=True)
args = parser.parse_args()

df = pd.DataFrame(columns=['Basename', 'sample_id', 'sample_well', 'sample_plate', 'sentrix_ID', 'sentrix_position', 'time', 'array_type', 'gender_donor', 'gender_recipient', 'eGFR_1month', 'eGFR_12month', 'eGFR_24month', 'DGF_status', 'DCD_status'])
for x in range(args.start, args.end): #Number of sheets within excel file 
	sheet_df =  pd.read_excel(args.excel, sheet_name=x, keep_default_na=True)
	sheet_df.columns = [col.strip() for col in sheet_df.columns] #Stripping whitespace in column names
	if x == 4: #This needed to be changed before creating the sample_id
		x = 12
	sample_id = 'SampleID_K' + str(x)
	for index in sheet_df.index: #iterate through columns
		if (sheet_df['Time'][index] == args.time1) or (sheet_df['Time'][index] == args.time2): #Only append to df if either time point is present 
			basename = str(sheet_df['Sentrix_ID'][index]) + '_' + str(sheet_df['Sentrix_Position'][index])
			try:
				row = sheet_df.loc[index, [sample_id, 'Sample_Well', 'Sample_Plate', 'Sentrix_ID', 'Sentrix_Position', 'Time', 'Array_type', 'gender_donor', 'gender_recipient', 'eGFR_1mo', 'eGFR_12mo', 'eGFR_24mo', 'DGF=1, non-DGF=0', 'DCD (yes/no)']]
				row = list(row)
				row.insert(0, basename)
				df.loc[len(df) + 1,] = row #adding to df
			except:
				row = sheet_df.loc[index, [sample_id, 'Sample_Well', 'Sample_Plate', 'Sentrix_ID', 'Sentrix_Position', 'Time', 'Array_type', 'gender_donor', 'gender_recipient', 'eGFR_1mo', 'eGFR_12mo', 'eGFR_24mo']]
				row = list(row)
				row.extend(['NaN', 'NaN'])
				row.insert(0, basename)
			finally:
				df.loc[len(df) + 1,] = row #adding to df
		else:
			pass

#print(arrayDict)

df = df[df['time'].notna()] #drop if time is empty

if args.month == '1':
	month = 'eGFR_1month'
	df.drop(['eGFR_12month', 'eGFR_24month', 'DGF_status', 'DCD_status'], axis=1, inplace=True)
elif args.month == '12':
	month = 'eGFR_12month'
	df.drop(['eGFR_1month', 'eGFR_24month', 'DGF_status', 'DCD_status'], axis=1, inplace=True)
elif args.month == '24':
	month = 'eGFR_24month'
	df.drop(['eGFR_1month', 'eGFR_12month', 'DGF_status', 'DCD_status'], axis=1, inplace=True)
else:
	print('Month not found')
	quit()

for index in df.index:
	if df[month][index] < 45:
		df[month][index] = 'Low'
	elif df[month][index] >= 45:
		df[month][index] = 'High'
	else:
		df[month][index] == 'NaN'

df.rename(columns={month:'eGFR'}, inplace=True)

os.chdir(args.working)
print(df)
df.to_csv(args.output)
