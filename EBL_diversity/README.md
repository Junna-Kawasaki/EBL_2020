# Assessment of genetic diversity of ancient bornaviral sequences   

## Discription  
To calculate genetic distance between ancient bornaviral sequences, please try the program: **EBL_diversity_batch.sh**.  

- This program provides an example to calculate genetic distance between ancient bornaviral N protein sequences.  
- Output file is **output_examples/genetic_distance.txt**.  

## Scripts  
1. EBL_genetic_distance.py  
- We used the most recent common ancestor (MRCA) of the EBLN orthologs as the ancestral bornaviral N gene to avoid overestimating the sequence diversity of ancient bornaviruses.  
- To provide a comparison standard for interpreting the ancient bornaviral genetic diversity, we calculated the genetic distance for classifying extant bornaviral species in our phylogenetic tree.  

Usage:  
python EBL_genetic_distance.py \  
	${tree_file} \  
	${age_file} \  
	${virus_name_file} \  
	${output_name}  

Parameters:  
- tree_file: phylogenetic tree file using EBLNs and extant bornaviral N proteins  
- age_file: file about EBL orthologous groups used to determine the MRCA node of each group  
- virus_name_file: virus name file to replace node name with virus species name  
- output_name: output file name  