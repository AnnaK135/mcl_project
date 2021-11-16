#!/bin/bash 

### For commenting out pieces of the code ###
[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

path_to_project="/mnt/e/MCL_project/chipseq_h3_ubf/"
path_to_ref="$path_to_project""ref/"
path_to_bwa_files="$path_to_project""bwa/"
path_to_rose_files="$path_to_project""rose/"

samples=(L1b L2b L3b Linput G1b G2b G3b Ginput)
lcl=(L1b L2b L3b)
granta=(G1b G2b G3b)

### Overlapping enhancer sets for GRANTA and control

# regular subtraction
bedtools subtract -a "$path_to_chipr"GRANTA_H3K27Ac_chipr_all.gff -b "$path_to_chipr"LCL_H3K27Ac_chipr_all.gff > granta_unique_H3K27Ac.gff
bedtools subtract -a "$path_to_chipr"LCL_H3K27Ac_chipr_all.gff -b "$path_to_chipr"GRANTA_H3K27Ac_chipr_all.gff > granta_unique_H3K27Ac.gff
# removing the whole feature if any overlap (-A)
bedtools subtract -A -a "$path_to_chipr"GRANTA_H3K27Ac_chipr_all.gff -b "$path_to_chipr"LCL_H3K27Ac_chipr_all.gff > granta_unique_H3K27Ac_stringent.gff
bedtools subtract -A -a "$path_to_chipr"LCL_H3K27Ac_chipr_all.gff -b "$path_to_chipr"GRANTA_H3K27Ac_chipr_all.gff > granta_unique_H3K27Ac_stringent.gff

### Overlapping super-enhancer sets for GRANTA and control 
