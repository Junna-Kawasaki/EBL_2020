#!/usr/bin/env python

import sys, re, os
import pandas as pd
import numpy as np
argvs = sys.argv

# 0. inputs and parameters
path = argvs[1]
gaps = argvs[2]
prop = argvs[3]

site_coverage = 1- float(gaps) # exclude sites where over ${gaps} of sequences were gaps
gap_thresh = 20 # Upper limit of consecutive gaps
len_block = 10 # Lower limit for consecutive alignment sites
prop_thresh = float(prop) #  remove sequences with less than ${prop} of the total alignment sites

# OPTION. specify sequences that you do not want to include in the coverage calculation.
try:
	exclude_f = argvs[4]
	exclude_l = []
	for line in open(exclude_f):
		line = line.strip()
		if ">" in line:
			exclude_l.append(line)
except:
	pass

# 1. read multiple alignement sequences
fasta_d = {}
table_d = {}
header_l = []
species_input = []
for line in open(path):
	line = line.strip()
	if ">" in line:
		header = line
		header_l.append(header)
		if header in exclude_l:
			fasta_d[header] = []
		else:
			table_d[header] = []
		species = line.split(";")[0]
		species = species.replace(">","")
		species_input.append(species)
	else:
		seq = line
		seq_l = [s for s in seq]
		if header in fasta_d.keys():
			fasta_d[header].extend(seq_l)
		else:
			table_d[header].extend(seq_l)


# 2. Calculate the gap proportion in each site
df = pd.DataFrame(table_d)
df_coverage = (df != "-")
df["coverage"] = df_coverage.sum(axis=1) / len(table_d)
df["seq_no"] = df_coverage.sum(axis=1)
df = df[df["coverage"] != 0.0]
df = df.reset_index()
df = df.T

# 3. Exclude sites where over ${gaps} of sequences were gaps
alignment_l = []
seq_no_l = df.loc["seq_no",:].values.tolist()
for site in range(0, df.shape[1] ):
	if df.at["seq_no", site] > site_coverage * max(seq_no_l):
		alignment_l.append(site)
		
# 4. Check continuity of alignment sites
l = []
continuous_l = []
for n in range(0,len(alignment_l)-1):
	in1 = n
	in2 = n+1
	if sorted(alignment_l)[in2] == sorted(alignment_l)[in1] +1:
		if in1 not in l:
			l.append(alignment_l[in1])
		l.append(alignment_l[in2])
	else:
		if len(l) > 0:
			continuous_l.append(l)
		l = []
	if in2 == len(alignment_l)-1:
		continuous_l.append(l)

##
aligned_site_l = []
continuous_l = [c for c in continuous_l if len(c) > 0]
if len(continuous_l) > 1:
	for n in range(0,len(continuous_l)-1):
		in1 = n
		in2 = n+1
		if n == 0:
			l = []
		if len(set(continuous_l[in1])) > len_block:
			gap = sorted(set(continuous_l[in2]))[0] - sorted(set(continuous_l[in1]))[-1]
			if gap < gap_thresh:
				if len(l) == 0:
					continuous_l[in1].extend(continuous_l[in2])
					l = continuous_l[in1]
				else:
					l.extend(continuous_l[in2])
			else:
				aligned_site_l.append(l)
				l = []
		if in2 == len(continuous_l)-1:
			aligned_site_l.append(l)
else:
	aligned_site_l = continuous_l[0]

##
aligned_site_l = [c for c in aligned_site_l if len(c) > 0]
length = 50
if len(aligned_site_l) == 0:
	for l in continuous_l:
		if len(set(l)) > length:
			length = len(set(l))
			max_l = l
	aligned_site_l.append(max_l)

##
length = 0
seq_inf = {}
if len(aligned_site_l) > 0:
	for n, l in enumerate(aligned_site_l):
		if len(l) > 1:
			seq_inf[n] = [sorted(set(l))[0], sorted(set(l))[-1], len(set(l))]
		if length < len(set(l)):
			length = len(set(l))
			max_n = n
elif len(aligned_site_l) == 1:
	length = len(set(aligned_site_l[0]))
	if length > 50:
		max_n = 0
	else:
		print ("There is No continuous alignment: " + str(length) + "aa")
		max_n = "NA"
else:
	print ("There is No continuous alignment")


## 5. output (log file and sequence file)
if length == 0:
	print ("No alignment with enough lengths")
elif length != 0 and max_n !="NA":
	aligned_site = aligned_site_l[max_n]
	
	df_aligned = df.loc[:, sorted(set(aligned_site))]
	virus_aligned_site = df_aligned.loc["index",:].values.tolist()
	df_aligned_prop = (df_aligned != "-")
	df_aligned["prop"] = df_aligned_prop.sum(axis=1) / df_aligned.shape[1]
	df_output = df_aligned[df_aligned["prop"] > prop_thresh]
	output_EBL_l = df_output.reset_index().values.tolist()
	if len(fasta_d) > 0:
		df_virus = pd.DataFrame(fasta_d)
		df_virus = df_virus.T
		df_virus = df_virus.loc[:, sorted(set(virus_aligned_site))]
		output_virus_l = df_virus.reset_index().values.tolist()
		output_l = output_virus_l + output_EBL_l
	else:
		output_l = output_EBL_l
	
	##
	o = "_".join(map(str,["", str(gaps), str(prop), "trimed.fasta"]))
	out = path.replace(".fasta", o)
	f = open(out,"w")
	for l in output_l: #output_rmdup: #output_sorted_l:
		if ">" in l[0] and l[0] not in exclude_l:
			gap_len = [s for s in l[1:-1] if s == "-"]
			if len(gap_len) == len(l[1:-1]):
				pass
			else:
				f.write(l[0] + "\n")
				f.write("".join(map(str,l[1:-1])) + "\n")
		elif ">" in l[0]:
			gap_len = [s for s in l[1:] if s == "-"]
			if len(gap_len) == len(l[1:]):
				pass
			else:
				f.write(l[0] + "\n")
				f.write("".join(map(str,l[1:])) + "\n")
	f.close()

