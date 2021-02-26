#!/bin/bash

input="input_examples/"
output="output_examples/repeatmasker/"
bed_f="hsEBLN1.bed"
fasta_f="NC_000010.11.fasta"

file=$(echo ${bed_f} | cut -d "." -f 1)
species="Homo sapiens"
flanking=10

mkdir ${output}

# 1. Extract nucleotide sequences 
## downstream
seqkit subseq \
	--bed  ${input}${bed_f} \
	${input}${fasta_f} \
	-d $((${flanking} * 1000)) -f \
	> ${output}${file}_${flanking}kb_fd_nucl.fasta

## upstream
seqkit subseq \
	--bed  ${input}${bed_f} \
	${input}${fasta_f} \
	-u $((${flanking} * 1000)) -f \
	> ${output}${file}_${flanking}kb_fu_nucl.fasta

## 
cat \
	${output}${file}_${flanking}kb_fd_nucl.fasta \
	${output}${file}_${flanking}kb_fu_nucl.fasta \
	> ${output}${file}_${flanking}kb_nucl.fasta 

# 2. Remove repetitive elements 
RepeatMasker \
	-q -xsmall -a -species ="${species}" \
	${output}${file}_${flanking}kb_nucl.fasta \
	-dir ${output}

python \
	fasta_remove_lower_case.py \
	${output}${file}_${flanking}kb_nucl.fasta.masked
	${output}${file}_${flanking}kb_nucl_RM.fasta