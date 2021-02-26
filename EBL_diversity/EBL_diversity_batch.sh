#!/bin/bash

input="input_examples/"
tree_file="query_out_EBLN_mafft_add_sed_cls_aa_filtering_100_rmdup_mafft_0.2_0.8_trimed.fasta.treefile"
age_file="query_out_EBL_mafft_add_sed_cls_aa_filtering_masked_RM_blastn_0.09_partition_sorted.txt"
output="output_examples/"

mkdir ${output}

# 1. Calculate genetic distances among nodes in a phylogenetic tree
python \
	EBL_genetic_distance.py \
	${input}${tree_file} \
	${input}${age_file} \
	${input}Bornavirus_species.txt \
	${output}genetic_distance.txt