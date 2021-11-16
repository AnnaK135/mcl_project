#!/bin/bash 

### For commenting out pieces of the code ###
[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

path_to_project="/mnt/e/MCL_project/chipseq_h3_ubf/"
path_to_ref="$path_to_project""ref/"
path_to_bwa_files="$path_to_project""bwa/"
path_to_rose_files="$path_to_project""rose/"
path_to_chipr="$path_to_project""chipr/"

samples=(L1b L2b L3b Linput G1b G2b G3b Ginput)
lcl=(L1b L2b L3b)
granta=(G1b G2b G3b)

### Overlapping enhancer sets for GRANTA and control
## for ChIPenrich  
# regular subtraction
bedtools subtract -a "$path_to_chipr"GRANTA_H3K27Ac_chipr_all.bed -b "$path_to_chipr"LCL_H3K27Ac_chipr_all.bed > "$path_to_chipr"granta_unique_H3K27Ac.bed
bedtools subtract -a "$path_to_chipr"LCL_H3K27Ac_chipr_all.bed -b "$path_to_chipr"GRANTA_H3K27Ac_chipr_all.bed > "$path_to_chipr"lcl_unique_H3K27Ac.bed

# another option to compare: removing the whole feature if any overlap (-A)
bedtools subtract -A -a "$path_to_chipr"GRANTA_H3K27Ac_chipr_all.bed -b "$path_to_chipr"LCL_H3K27Ac_chipr_all.bed > "$path_to_chipr"granta_unique_H3K27Ac.bed
bedtools subtract -A -a "$path_to_chipr"LCL_H3K27Ac_chipr_all.bed -b "$path_to_chipr"GRANTA_H3K27Ac_chipr_all.bed > "$path_to_chipr"lcl_unique_H3K27Ac.bed

# filtering out the small regions < 10 bp
awk '($3-$2) >= 10' "$path_to_chipr"granta_unique_H3K27Ac.bed > "$path_to_chipr"granta_unique_H3K27Ac_filtered10.bed
awk '($3-$2) >= 10' "$path_to_chipr"lcl_unique_H3K27Ac.bed > "$path_to_chipr"lcl_unique_H3K27Ac_filtered10.bed

# selecting chr11 and chr14
awk '($1) == "chr11"' "$path_to_chipr"granta_unique_H3K27Ac_filtered10.bed > "$path_to_chipr"granta_unique_H3K27Ac_filtered10_chr11.bed
awk '($1) == "chr14"' "$path_to_chipr"granta_unique_H3K27Ac_filtered10.bed > "$path_to_chipr"granta_unique_H3K27Ac_filtered10_chr14.bed
awk '($1) == "chr11"' "$path_to_chipr"lcl_unique_H3K27Ac_filtered10.bed > "$path_to_chipr"lcl_unique_H3K27Ac_filtered10_chr11.bed
awk '($1) == "chr14"' "$path_to_chipr"lcl_unique_H3K27Ac_filtered10.bed > "$path_to_chipr"lcl_unique_H3K27Ac_filtered10_chr14.bed

 
