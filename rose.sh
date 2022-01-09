#!/bin/bash
### For commenting out pieces of the code ###
[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

### Introducing variables for paths so we can change them later and reuse the script ###
path_to_project="/mnt/e/MCL_project/chipseq_h3_ubf/"
path_to_ref="/mnt/e/MCL_project/ref/"
path_to_tools="/mnt/e/tools/"
path_to_bwa_files="$path_to_project""bwa/"
path_to_rose_files="$path_to_project""rose/"

samples=(L1b L2b L3b Linput G1b G2b G3b Ginput)
lcl=(L1b L2b L3b)
granta=(G1b G2b G3b)

BEGINCOMMENT
### Indexing the BAM files
for i in ${samples[*]};do
    samtools index -b -@ 20 "$path_to_bwa_files"${i}_sorted_marked_filtered.bam &
done

### Calling superenhancers with ROSE ### 
cd "$path_to_tools"rose

for i in ${lcl[*]};do
    python ROSE_main.py -g hg38 -i "$path_to_rose_files"LCL_H3K27Ac_chipr_all.gff \
                                -r "$path_to_bwa_files"${i}_sorted_marked_filtered.bam \
                                -o "$path_to_rose_files"rose_output_${i} -t 2500 \
                                -c "$path_to_bwa_files"Linput_sorted_marked_duplicates.bam &
done

cd "$path_to_tools"rose
for i in ${granta[*]};do
    python ROSE_main.py -g hg38 -i "$path_to_rose_files"GRANTA_H3K27Ac_chipr_all.gff \
                                -r "$path_to_bwa_files"${i}_sorted_marked_filtered.bam \
                                -o "$path_to_rose_files"rose_output_${i} -t 2500 \
                                -c "$path_to_bwa_files"Ginput_sorted_marked_duplicates.bam &
done
ENDCOMMENT

chipr -i "$path_to_rose_files"rose_output_G1b/GRANTA_H3K27Ac_chipr_all_Gateway_SuperEnhancers.bed \
         "$path_to_rose_files"rose_output_G2b/GRANTA_H3K27Ac_chipr_all_Gateway_SuperEnhancers.bed \
         "$path_to_rose_files"rose_output_G3b/GRANTA_H3K27Ac_chipr_all_Gateway_SuperEnhancers.bed \
          -m 2 -o "$path_to_rose_files"GRANTA_SEs_chipr 