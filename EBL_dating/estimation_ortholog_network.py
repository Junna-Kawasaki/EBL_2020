#!/usr/bin/env python

### shirokaneで年代推定する

import sys, re
import os
import time
import networkx as nx
import community
import math

import numpy as np
import pandas as pd
pd.set_option("display.max_colwidth", 20)
pd.set_option("display.max_columns", 20)
pd.set_option("display.max_rows", 300)

# 0. input and parameters
argvs = sys.argv
blastn_f = argvs[1] 
blastn_cov = float(argvs[2]) 
out_path = argvs[3]

# 1. read blastn file
def blastn_reader(blastn_f):
	### read aligned sequence longer than 100 bp
	res_d = {}
	no_hit_n = []
	blastn_name_l = []
	for line in open(blastn_f):
		if '# Query:' in line:
			line = line.strip().split(" ")
			name = line[2]
			name = ";".join(name.split(";")[:-1])
			blastn_name_l.append(name)
		if '# 0 hits' in line:
			no_hit_n.append(name)
		if '#' not in line:
			line = line.strip().split('\t')
			query,subject,a_length,q_length,s_length,mismatch,gap,q_start,q_end,s_start,s_end,evalue = line[:12]
			query = re.sub(r'(.+)gi\|[0-9]+\|ref\|(.+)\|(.+)',r'\1\2\3',query)
			subject = re.sub(r'(.+)gi\|[0-9]+\|ref\|(.+)\|(.+)',r'\1\2\3',subject)
			evalue = float(evalue)
			a_length,q_length,s_length,mismatch,gap,q_start,q_end,s_start,s_end = [int(c) for c in [a_length,q_length,s_length,mismatch,gap,q_start,q_end,s_start,s_end]]
			query,subject = [str(c) for c in [query,subject]]
			if s_start < s_end:
				strand = '+'
			else:
				strand = '-'
			line.append(strand)
			if query not in res_d:
				res_d[query] = {}
			if subject not in res_d[query]:
				res_d[query][subject] = []
			if (a_length > 100):
				l = [query,subject,a_length,q_length,s_length,q_start,q_end,s_start,s_end,evalue,strand]
				res_d[query][subject].append(l)
	
	### summary alignment lengths
	blastn_d = {}
	alen_d = {}
	for name in blastn_name_l:
		alen_d[name] = {}
		for n in blastn_name_l:
			alen_d[name][n] = 1
	
	for query_name in res_d.keys():
		if query_name not in blastn_d.keys():
			blastn_d[query_name] = {}
		for subject_name in res_d[query_name].keys():
			if 'ups' in subject:
				flanking = 'ups'
			if 'dws' in subject:
				flanking = 'dws'
			if (len(res_d[query_name][subject_name]) == 1):
				query,subject,a_length,q_length,s_length,q_start,q_end,s_start,s_end,evalue,strand = res_d[query_name][subject_name][0]
				if (a_length > 10000 * blastn_cov):
					blastn_d[query_name][subject_name] = [q_length,a_length,flanking]
			else:
				index_name = ['query','subject','a_length','q_length','s_length','q_start','q_end','s_start','s_end','evalue','strand']
				res_l = res_d[query_name][subject_name] # res_d[query_name][subject_name]
				col_name = [str(c) for c in range(0,len(res_d[query_name][subject_name]))]
				df = pd.DataFrame(res_l, columns = index_name)
				#df =df.T
				df = df.sort_values(by=['a_length'])
				pd.set_option('display.width', 300)
				pd.set_option('display.max_columns',100)
				#print (df)
				qname = ";".join(query_name.split(";")[:-1])
				sname = ";".join(subject_name.split(";")[:-1])
				alen_d[qname][sname] = math.log10(a_length)
				if (df['a_length'].sum() > 10000 * blastn_cov):
					a_length = df['a_length'].sum() 
					blastn_d[query_name][subject_name] = [q_length,a_length,flanking]
	
	### edge between subject and query
	name_l_l = []
	blastn_l = []
	blastn_d2 = {}
	for query in blastn_d.keys():
		query_name = query
		#query_name = re.sub(r'gi\|[0-9]+\|ref\|([A-Za-z0-9]\.[0-9])\|',r'\1',query_name)
		query_name =  re.sub(r'(.+);dws',r'\1',query_name)
		query_name =  re.sub(r'(.+);ups',r'\1',query_name)
		name_l_l.append(query_name)
		if query_name not in blastn_d2.keys():
			blastn_d2[query_name] = []
		for subject in blastn_d[query].keys():
			subject_name = subject
			#subject_name = re.sub(r'gi\|[0-9]+\|ref\|([A-Za-z0-9]\.[0-9]+)\|',r'\1',subject_name)
			subject_name = re.sub(r'(.+);dws',r'\1',subject_name)
			subject_name = re.sub(r'(.+);ups',r'\1',subject_name)
			l = (query_name,subject_name)
			#print (l)
			blastn_l.append(tuple(sorted(l)))
			if subject_name not in blastn_d2[query_name]:
				blastn_d2[query_name].append(subject_name)
	
	return (alen_d, blastn_l)
	
# 2. Network analysis
def network_blast(ortholog_l):
	ortho_l = [list(i) for i in ortholog_l]
	sum_ortho_l = sum(ortho_l,[])
	len_ortho_l = len(set(sum_ortho_l))
	G = nx.Graph()
	G.add_edges_from(ortholog_l)
	subgroup_l = sorted([set(x) for x in nx.connected_components(G)])
	partition_d = community.best_partition(G)
	
	return (partition_d, G)

### 3. Output
def output_result(partition_d, G):
	partition_d_d = {}
	partition_l = []
	for c in partition_d.keys():
		n = partition_d[c]
		if n not in partition_d_d:
			partition_d_d[n] = []
			partition_d_d[n].append(c)
			partition_l.append(c)
		else:
			partition_d_d[n].append(c)
			partition_l.append(c)
	
	##
	for name in alen_d.keys():
		if name not in partition_d.keys():
			n += 1
			partition_d_d[n] = [name]
	
	##
	if os.path.exists(out_path):
		os.remove(out_path)
	
	## 
	group_l = []
	f = open(out_path, 'w')
	n = -1
	labels = {}
	labels_d ={}
	lis = []
	for c in sorted(partition_d_d.keys()):
		for i in sorted(partition_d_d[c]):
			n = n+1
			f.write ("\t".join(map(str, [n, "Group" + str(c), i])) + "\n")
			group_l.append(i)
			if i in G.nodes():
				labels[i] = n
			if c not in labels_d.values():
				labels_d[i] = c
	
	f.close()

# --------------------------------------------------
###  ---- Main -----
# 1. Read blastn file
## Run time
start = time.time()

alen_d, blastn_l = blastn_reader(blastn_f)
##
print (" --- finish read blastn ---")
elapsed_time = time.time() - start
print ("# elasped_time: " + str(elapsed_time))
# --------------------------------------------------
# 2. Network analysis
## Run time
start = time.time()

ortholog_l =  set(blastn_l) 
partition_d, G = network_blast(ortholog_l)

##
print (" --- finish network analysis ---")
elapsed_time = time.time() - start
print ("# elasped_time: " + str(elapsed_time))
# --------------------------------------------------
# 3. Output
## Run time
start = time.time()

output_result(partition_d, G)
##
print (" --- finish network analysis ---")
elapsed_time = time.time() - start
print ("# elasped_time: " + str(elapsed_time))
# --------------------------------------------------