#!/bin/bash

input="input_examples/"
output="output_examples/network/"
file="EBLN1-3_10kb_RM"

mkdir ${output}

# 1. Network analysis
criteria=0.09
python estimation_ortholog_network.py \
	${input}${file}_blastn \
	${criteria} \
	${output}${file}_network.txt