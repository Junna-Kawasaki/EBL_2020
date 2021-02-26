# Identification of EBLs in vertebrate genomic data  

## Description  
To reconstruct EBL sequences from tBLASTn search results, please try the program: **EBL_reconstruction_batch.sh**.  

- This program provides an example to obtain amino acid sequences of EBLNs in the human genome.  
- Output file is **output_examples/tBLASTn_human_EBLN_mafft_add_cls_aa.fasta**.  

## Scripts  
1. EBL_reconstruction.py  
Because most EBLs were detected as fragmented sequences due to mutations occurring after integration, we reconstructed EBL sequences by concatenating sequences using this script.  
