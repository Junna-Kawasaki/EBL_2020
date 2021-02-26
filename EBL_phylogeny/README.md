# Phylogenetic analysis  

## Description  
To costruct phylogenetic tree using extant viruses and EBLs, please try the program: **EBL_phylogenetic_analysis_batch.sh**  

- This program provides an example to construct phylogenetic tree using N proteins of extant viruses and EBLNs.  
- Output file is **output_examples/query_out_EBLN_mafft_add_sed_cls_aa_filtering_100_rmdup_mafft_0.3_0.7_trimed.fasta.treefile**.  

## Scripts  
1. alignment_seq_num.py  
By using this script, we can trim a multiple sequence alignment by excluding sites where over ${gaps}% of sequences were gaps, and subsequently removing sequences with less than ${site_coverage}% of the total alignment sites.  

Usage:  
python alignment_seq_num.py \  
	${multiple_sequence_alingment} \  
	${gaps} ${site_coverage} \  
	${exclude_sequences}  

Parameters:  
- multiple_sequence_alignment: input fasta file that you want to trim  
- gaps: proportion of gaps  
e.g. if you set this as 0.3, this script removes sites where over 30% of sequences were gaps.  
- site_coverage: proportion of aligned sites.  
e.g. if you set this as 0.7, this script removes sequences with less than 70% of the total alignment sites.  
- exclude_sequences  
This is an option. You can specify sequences that you do not want to include in the coverage calculation.    