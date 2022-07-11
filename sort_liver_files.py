#! /usr/bin/env python3

#Example input: ./sort_liver_files.py -b <base_dir> -o <out_dir>

import pandas as pd
import argparse as ap
import os
import shutil 

parser = ap.ArgumentParser()
parser = ap.ArgumentParser(description='Sort raw liver files')
parser.add_argument("-b", "--base_dir", help="base directory", required=True)
parser.add_argument("-o", "--out_dir", help="output directory", required=True)
args = parser.parse_args()

#Move to working directory and sort files into respective folders 
os.chdir(args.base_dir)

#Append to dir_list
dir_list = []
for directory in os.listdir(args.base_dir):
	directory = os.path.join(args.base_dir, directory)
	dir_list.append(directory)

#Iterate through directories
file_list = []
for direct in dir_list:
	for data in os.listdir(direct):
		if ".idat" in data:
			data = os.path.join(direct, data)
			file_list.append(data)


#print(file_list)
os.chdir(args.out_dir)

if os.path.isdir(args.out_dir + 'control/') and os.path.isdir(args.out_dir + '450k_array/') and os.path.isdir(args.out_dir + 'EPIC_array/'):
	print('Skip making directories')
else:
	os.makedirs('control')
	os.makedirs('450k_array')
	os.makedirs('EPIC_array')

#Copying to output directory
for file in file_list:
	src = file
	if "450K" in src:
		extension = file.split("/")
		extension = extension[-1]
		dst = args.out_dir + '450k_array/' + extension
		shutil.copyfile(src, dst)
	elif "EPIC" in src:
		extension = file.split("/")
		extension = extension[-1]
		dst = args.out_dir + 'EPIC_array/' + extension
		shutil.copyfile(src, dst)
	else:
		extension = file.split("/")
		extension = extension[-1]
		dst = args.out_dir + 'control/' + extension
		shutil.copyfile(src, dst)
	



