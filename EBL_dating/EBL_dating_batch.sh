#!/bin/bash

input="input_examples/"
output="output_examples/"

mkdir ${output}

# 1. Compare flanking sequences of EBLs
## remove repetitive elements by repeatmasker
bash estimation_age_repeatmasker.sh

## compare similarity among flanking sequences by blastn
bash estimation_age_blastn.sh

# 2. Network analysis
bash estimation_age_network.sh
