#!/usr/bin/env python

import sys, re
argvs = sys.argv

fasta = open(argvs[1])
out_name = argvs[2]

seq_d = {}
for line in fasta:
	line = line.strip()
	if '>' in line:
		name=line.replace('>','')
		seq_d[name] = ''
	else:
		seq = line
		seq = re.sub(r'a|t|c|g','',seq)
		seq_d[name] = seq_d[name] + seq

out_f = out_name + '_RM.fasta'
f = open(out_f, 'w')
for name in seq_d.keys():
	seq = seq_d[name]
	f.write('>' + name + '\n')
	f.write(seq + '\n')
f.close()