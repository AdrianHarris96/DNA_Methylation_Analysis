#! /usr/bin/env python3

import pandas as pd
import os

df = pd.read_csv('/Users/adrianharris/Documents/dna_methylation_analysis/kidney_sample_sheet.csv')

sample_dictionary = {}
unknown_dictionary ={}
for index in df.index:
	sample = df['sample_id'][index]
	sample = sample.strip()
	if sample[-3] == 'K' or sample[-3] == 'k' or sample[-2] == 'K' or sample[-2] == 'k':
		if sample[:-1] not in sample_dictionary.keys():
			sample_dictionary[sample[:-1]] = 1
		else:
			sample_dictionary[sample[:-1]] += 1
	else:
		if sample not in unknown_dictionary.keys():
			unknown_dictionary[sample] = 1
		else:
			unknown_dictionary[sample] += 1

sample_dictionary = dict(sorted(sample_dictionary.items()))

final_dictionary = {}
for key in sample_dictionary.keys():
    if sample_dictionary[key] == 2:
        final_dictionary[key] = 2

print(final_dictionary)
print(len(final_dictionary))

paired_df = pd.DataFrame(columns = df.columns)

for index in df.index:
	sample = df['sample_id'][index]
	if sample[:-1] in final_dictionary.keys():
		paired_df.loc[len(paired_df)] = list(df.loc[index,])

paired_df = paired_df.iloc[: , 1:]
print(paired_df)
paired_df.to_csv('/Users/adrianharris/Desktop/paired_kidney_sample_sheet.csv', index = False)