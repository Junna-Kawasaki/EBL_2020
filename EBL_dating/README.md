# Dating analysis for the integrated age of EBLs  

## Discription  
To estimate EBL integration ages by determining orthologous loci, please try the program: **EBL_dating_batch.sh**.  

- This program provides an example to cluster EBLNs into orthologous groups: hsEBLN-1 to hsEBLN-3.  
- Output file is **output_examples/network/EBLN1-3_10kb_RM_network.txt**.  

## Scripts  
1. estimation_age_repeatmasker.sh  
This program provides an example to remove reptitive elements from the downstream and upstream sequences of hsEBLN-1.  
** Please causion** 
To try this program, please prepare fasta file of "NC_000010.11.fasta" in directory: "input_examples/.  

2. estimation_age_blastn.sh  
This program provides an example to compare sequence similarities of the EBLN flanking regions by BLASTN.  

- fasta_remove_lower_case.py  
The lower cases in a fasta file can be removed by this script.  
This program requires a fasta file from which you want to remove lower cases.  

4. estimation_age_network.sh  
This program provides an example to cluster orthologous EBLN loci based on the detection of community structure in the sequence similarity network.  

5. estimation_ortholog_network.py  
This script performes making sequence similarity network and extracting community structures.

Usage: 
python estimation_ortholog_network.py  \  
	${blastn_file} \  
	${criteria} \  
	${output_file}  
	
Parameters:  
- blastn_file: input file obtained from all-against-all sequence comparison by BLASTN  
- criteria: criteria to construct a sequence network  
e.g. If we set this value as 0.9, network nodes were connected when their sequence alignment coverage was over 9.0% of the entire sequence length.  
- output_file: name of output file
