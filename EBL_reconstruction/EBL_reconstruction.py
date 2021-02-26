#!/usr/bin/env python

#### python 3
import sys, re
import os
import numpy as np
import pandas as pd
pd.set_option("display.max_colwidth", 20)
pd.set_option("display.max_columns", 500)
pd.set_option("display.max_rows", 300)
pd.options.display.precision = 2
pd.options.display.float_format = '{:.2e}'.format

argvs = sys.argv

# 0. input files
mafft_f = argvs[1] 
mapout_f = argvs[2] 

# defined criteria to concate bornavirus-like sequences
domain = re.sub(r".+EBL(.)_.+",r"\1",mafft_f)
## [genomic distance, distance of aligned positions to bornaviral proteins]
threshold_d = {"N": [1000, 300], "P": [1000, 300], "M":[1000, 300], "G":[1000,300], "L":[2000, 650]}
g_thresh, q_thresh = threshold_d[domain]

# 1. function of read fasta file
def fasta_reader(path,dic):
	for line in open(path):
		line = line.strip()
		if ">" in line:
			header = line
			if header not in dic:
				dic[header] = []
		else:
			seq = [s for s in line]
			dic[header].extend(seq)
	return (dic)

# 2. read mafft alignment file (fasta)
mafft_d = {}
mafft_d = fasta_reader(mafft_f, mafft_d)

# 3. read information file on alignment positions to bornavirus proteins
mapout_d = {}	
for line in open(mapout_f):
	if ">" in line:
		name = line.strip()
		residue_l = []
		original_l = []
		ref_l = []
	if "#" in line:
		reference = line
	if (">" not in line) and ("#" not in line):
		info = line.strip().split(", ")
		residue,original,ref = info[:3]
		residue_l.append(residue)
		original_l.append(original)
		ref_l.append(ref)
	mapout_d[name] = [residue_l, original_l, ref_l]

# 4. read tblastn information
info_pos_d = {}
info_neg_d = {}
for key in mapout_d.keys():
	header = key
	key = key.replace(">","")
	species, chrom, start,end,query,evalue,qstart,qend,strand = key.split(";")[:9]
	original_l, ref_l = mapout_d[header][1:]
	ref_l = [int(r) for r in ref_l if r != "-"]
	ref_start = min(ref_l)
	ref_end = max(ref_l)
	if strand == "+":
		info_pos_d[header] = [species, chrom, start,end,query,evalue,ref_start,ref_end,strand]
	else:
		info_neg_d[header] = [species, chrom, end,start,query,evalue,ref_start,ref_end,strand]

## 5. MERGE: compared genomic distance and aligned position with criteria
dis_pos_mrg_l = []
dis_neg_mrg_l = []

# positive strand
df_pos = pd.DataFrame(info_pos_d, index= ["species", "chrom", "start","end","query","evalue","ref_start","ref_end","strand"])
df_pos = df_pos.T
df_pos = df_pos.astype({"start": int, "end":int, "ref_start":int, "ref_end":int})
df_pos = df_pos.sort_values(["chrom","start"], ascending=[True,True])

mrg_pos_l = []
for i in range(0, df_pos.shape[0]-1):
	key1 = df_pos.index.tolist()[i]
	key2 = df_pos.index.tolist()[i+1]
	if df_pos.loc[key1,"chrom"] == df_pos.loc[key2,"chrom"]:
		g_dis = df_pos.loc[key2,"start"] - df_pos.loc[key1,"end"]
		q_dis = df_pos.loc[key2,"ref_start"] - df_pos.loc[key1,"ref_end"]
		if g_dis < 1: #and q_dis < 1:
			l = [i, i+1]
			mrg_pos_l.append(l)
			dis_pos_mrg_l.append([g_dis,q_dis])
		else:
			mrg_pos_l.append([i])
			mrg_pos_l.append([i+1])
			dis_pos_mrg_l.append([g_dis,q_dis])
	else:
		mrg_pos_l.append([i])
		mrg_pos_l.append([i+1])
		dis_pos_mrg_l.append(["-","-"])

dis_pos_mrg_l.append(["***","***"])

## negative strand
df_neg = pd.DataFrame(info_neg_d, index= ["species", "chrom", "start","end","query","evalue","ref_start","ref_end","strand"])
df_neg = df_neg.T
df_neg = df_neg.astype({"start": int, "end":int, "ref_start":int, "ref_end":int})
df_neg = df_neg.sort_values(["chrom","end"], ascending=[True,False])

mrg_neg_l = []
for i in range(0, df_neg.shape[0]-1):
	key1 = df_neg.index.tolist()[i]
	key2 = df_neg.index.tolist()[i+1]
	if df_neg.loc[key1,"chrom"] == df_neg.loc[key2,"chrom"]:
		g_dis = df_neg.loc[key1,"end"] - df_neg.loc[key2,"start"]
		q_dis = df_neg.loc[key2,"ref_start"] - df_neg.loc[key1,"ref_end"]
		if g_dis < 1: #g_dis < 1 and q_dis < 1:
			l = [i, i+1]
			mrg_neg_l.append(l)
			dis_neg_mrg_l.append([g_dis,q_dis])
		else:
			mrg_neg_l.append([i])
			mrg_neg_l.append([i+1])
			dis_neg_mrg_l.append([g_dis,q_dis])
	else:
		dis_neg_mrg_l.append(["-","-"])
		mrg_neg_l.append([i])
		mrg_neg_l.append([i+1])

dis_neg_mrg_l.append(["***","***"])

## 6. MERGE sequences
def list_concat(list_l):
	concat_l = [set(list_l[0])]
	for s in list_l[1:]:
		last_s = concat_l[-1]
		if len(set(last_s) & set(s)) > 0:
			concat_l[-1] = set(last_s) | set(s)
		else:
			concat_l.append(set(s))
	return (concat_l)

mrg_pos_l = list_concat(mrg_pos_l)
mrg_neg_l = list_concat(mrg_neg_l)

## 7. CONCAT: compared genomic distance and aligned position with criteria
dis_pos_cls_l = []
dis_neg_cls_l = []

## positive strand
cls_pos_l = []
for m in range(0, len(mrg_pos_l)-1):
	df1 = df_pos.iloc[list(mrg_pos_l[m]),:]
	df2 = df_pos.iloc[list(mrg_pos_l[m+1]),:]
	g_dis = "-"
	q_dis = "-"
	if df1["chrom"].values.tolist()[0] == df1["chrom"].values.tolist()[0]:
		g_dis = min(df2["start"].values.tolist()) - max(df1["end"].values.tolist())
		q_dis = min(df2["ref_start"].values.tolist()) - max(df1["ref_end"].values.tolist())
		if -g_thresh < g_dis < g_thresh and -g_thresh<q_dis <g_thresh:
			l = list(mrg_pos_l[m]) + list(mrg_pos_l[m+1])
			cls_pos_l.append(l)
			dis_pos_cls_l.append([g_dis,q_dis])
		else:
			cls_pos_l.append(list(mrg_pos_l[m]))
			cls_pos_l.append(list(mrg_pos_l[m+1]))
			dis_pos_cls_l.append([g_dis,q_dis])
	else:
		dis_pos_cls_l.append(["-","-"])
		cls_pos_l.append(list(mrg_pos_l[m]))
		cls_pos_l.append(list(mrg_pos_l[m+1]))

dis_pos_cls_l.append(["***","***"])

## negative strand
cls_neg_l = []
for m in range(0, len(mrg_neg_l)-1):
	df1 = df_neg.iloc[list(mrg_neg_l[m]),:]
	df2 = df_neg.iloc[list(mrg_neg_l[m+1]),:]
	g_dis = "-"
	q_dis = "-"
	if df1["chrom"].values.tolist()[0] == df1["chrom"].values.tolist()[0]:
		g_dis = min(df1["end"].values.tolist()) - max(df2["start"].values.tolist())
		q_dis = min(df2["ref_start"].values.tolist()) - max(df1["ref_end"].values.tolist())
		if -g_thresh < g_dis < g_thresh and -g_thresh<q_dis <g_thresh: 
			l = list(mrg_neg_l[m]) + list(mrg_neg_l[m+1])
			cls_neg_l.append(l)
			dis_neg_cls_l.append([g_dis,q_dis])
		else:
			cls_neg_l.append(list(mrg_neg_l[m]))
			cls_neg_l.append(list(mrg_neg_l[m+1]))
			dis_neg_cls_l.append([g_dis,q_dis])
	else:
		dis_neg_cls_l.append(["-","-"])
		cls_neg_l.append(list(mrg_neg_l[m]))
		cls_neg_l.append(list(mrg_neg_l[m+1]))

dis_neg_cls_l.append(["***","***"])

cls_pos_l = list_concat(cls_pos_l)
cls_neg_l = list_concat(cls_neg_l)

### 8. sort bornavirus-like sequences according to E-value in tBLASTn search
df_mafft = pd.DataFrame(mafft_d)
df_mafft = df_mafft.T

## 
def output_seq(list_l, df):
	aa_d = {}
	except_l = []
	for l in list_l:
		df_s = df.iloc[list(l),:]
		df_s = df_s.astype({"start": int, "end":int, "ref_start":int, "ref_end":int, "evalue": float})
		df_s = df_s.sort_values("evalue", ascending = True)
		species = df_s["species"].values.tolist()[0]
		chrom = df_s["chrom"].values.tolist()[0]
		strand = df_s["strand"].values.tolist()[0]
		evalue = min(df_s["evalue"].values.tolist())
		evalue = f'{evalue:.2e}'
		query = df_s["query"].values.tolist()[0]
		ref_start = min(df_s["ref_start"].values.tolist())
		ref_end = max(df_s["ref_end"].values.tolist())
		if strand == "+":
			start = min(df_s["start"].values.tolist())
			end = max(df_s["end"].values.tolist())
		else:
			 start = min(df_s["end"].values.tolist())
			 end = max(df_s["start"].values.tolist())
		chrom = re.sub(r"gi\|[0-9]*\|ref\|([A-Z]*_[0-9]*.[0-9]*)\|",r"\1",chrom)
		header_l = [species, chrom, start, end, query, evalue, ref_start, ref_end, strand]
		header = ";".join(map(str, header_l))
		
		key_l = df_s.index.tolist()
		aa_l = ["-"] * df_mafft.shape[1]
		for k in key_l:
			for i in range(0, df_mafft.shape[1]):
				base = df_mafft.loc[k,i]
				if base != "-" and aa_l[i] == "-":
					aa_l[i] = base
				else:
					pass
					
		aa_d[header] = aa_l	
	
	return (aa_d, except_l)

##
mrg_pos_aa, except_pos_mrg = output_seq(mrg_pos_l, df_pos)
mrg_neg_aa, except_neg_mrg = output_seq(mrg_neg_l, df_neg)
cls_pos_aa,  except_pos_cls = output_seq(cls_pos_l, df_pos)
cls_neg_aa,  except_neg_cls = output_seq(cls_neg_l, df_neg)


### 9. output sequence file
def output_fasta(path, dic):
	f = open(path, "a")
	for key,val in dic.items():
		f.write(">" + key + "\n")
		val = [s for s in val if s != "-"]
		val = [s for s in val if s != "---"]
		f.write("".join(map(str,val)) + "\n")
	f.close()

cls_aa_out = mafft_f.replace(".fasta", "_cls_aa.fasta")
mrg_aa_out = out = mafft_f.replace(".fasta", "_mrg_aa.fasta")

##
for path in  [cls_aa_out, mrg_aa_out]:
	if os.path.exists(path):
			os.remove(path)

##
for path in [cls_aa_out, mrg_aa_out]:
		if "cls" in path and "aa" in path:
			output_fasta(path, cls_pos_aa)
			output_fasta(path, cls_neg_aa)
		elif "mrg" in path and "aa" in path:
			output_fasta(path, mrg_pos_aa)
			output_fasta(path, mrg_neg_aa)	
