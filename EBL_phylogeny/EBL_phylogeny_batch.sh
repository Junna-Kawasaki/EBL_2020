#!/bin/bash

input="input_examples/"
file="query_out_EBLN_mafft_add_sed_cls_aa_filtering"
virus_f="Bornaviridae_Nyami_protein_N_0.9"
output="output_examples/"
domain="N"
len=100

mkdir ${output}

# 1. Select EBL sequences 
## based on the length (longer than 100 aa)
seqkit seq -m ${len} \
	${input}${file}.fasta.gz \
	> ${output}${file}_${len}.fasta

## remove sequence redundancy
seqkit rmdup -s ${output}${file}_${len}.fasta \
	-o ${output}${file}_${len}_rmdup.fasta

# 2. make multiple sequence alignment	
## make alignment using virus proteins
mafft --auto ${input}${virus_f}.fasta  \
  > ${output}${virus_f}_mafft.fasta

## make alignment using virus proteins and EBLs
mafft --keeplength \
	--addfragments ${output}${file}_${len}_rmdup.fasta \
	--auto ${output}${virus_f}_mafft.fasta \
	> ${output}${file}_${len}_rmdup_mafft.fasta
	
# 3. Select aligned positions
gaps=0.3
prop=0.7
python alignment_seq_num.py \
	${output}${file}_${len}_rmdup_mafft.fasta \
	${gaps} \
	${prop} \
	${output}${virus_f}_mafft.fasta

# 4. Make phylogenetic tree
rm ${output}${file}_${len}_rmdup_mafft_${gaps}_${prop}_trimed.fasta.*

iqtree \
	-s ${output}${file}_${len}_rmdup_mafft_${gaps}_${prop}_trimed.fasta \
	-m TESTNEW \
	-bb 1000 \
	-alrt 1000 \
	-redo \
	-nt AUTO