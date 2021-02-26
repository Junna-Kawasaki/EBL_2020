#!/bin/bash

input="input_examples/"
output="output_examples/blastn/"
file="EBLN1-3_10kb_RM"

mkdir ${output}

# 1. Compare sequence similarities of EBLN flanking sequences
# Make database
makeblastdb \
	-in ${input}${file}.fasta \
    -dbtype nucl \
    -out ${output}${file}_DB

# BLASTN
blastn \
	-db ${output}${file}_DB \
	-query ${input}${file}.fasta \
	-out ${input}${file}_blastn \
	-evalue 0.01 \
	-outfmt "7 qseqid sseqid length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore pident"