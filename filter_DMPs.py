#!/usr/bin/env python3

#Example input: ./filter_DMPs.py -d /Users/adrianharris/Desktop/e1month/

import pandas as pd
import argparse as ap
import os

parser = ap.ArgumentParser()
parser = ap.ArgumentParser(description='filter DMPs.csv to DMPs_sig.csv')
parser.add_argument("-d", "--directory", help="path to directory of DMPs.csv", required=True)
args = parser.parse_args()

directory = args.directory
f = []
for filename in os.listdir(directory):
    file = os.path.join(directory, filename)
    f.append(file)

os.chdir(args.directory)

for file in f:
	print(f"Running file {file}")
	DMPs_df = pd.read_csv(file)
	DMPs_df = DMPs_df.loc[DMPs_df['adj.P.Val'] < 0.05]
	DMPs_df.sort_values(['adj.P.Val'], ascending=True, inplace=True)
	ouptut_name = file[:-4]
	ouptut_name += "_sig.csv"
	DMPs_df.to_csv(ouptut_name)
	