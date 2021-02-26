#!/bin/bash

input="input_examples/"
file="tBLASTn_human"
virus_f="input_examples/Bornaviridae_protein_sed_80"
output="output_examples/"
domain="N"

mkdir ${output}

## make sequence alignment by mafft using extant bornavirus protein sequences and bornavirus-like sequences detected in tBLASTn screening

mafft \
	--auto ${virus_f}_${domain}.fasta \
	>  ${virus_f}_${domain}_mafft.fasta

mafft --keeplength --mapout \
	--addfragments ${input}${file}_${domain}.fasta \
	--auto ${virus_f}_${domain}_mafft.fasta \
	> ${output}${file}_${domain}_mafft_add.fasta
  
mv ${input}${file}_${domain}.fasta.map ${output}

cat ${output}${file}_${domain}_mafft_add.fasta | \
  seqkit grep -v -n -r -p "virus" \
  > ${output}${file}_EBL${domain}_mafft_add.fasta

## reconstruct EBL sequence based on the genomic locations and alignment position to bornavirus proteins

python \
	EBL_reconstruction.py \
	${output}${file}_EBL${domain}_mafft_add.fasta \
	${output}${file}_${domain}.fasta.map 