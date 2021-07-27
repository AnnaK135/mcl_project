#!/bin/bash

### For commenting out pieces of the code ###
[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

### Introducing variables for paths so we can change them later and reuse the script ###
path_to_project="/mnt/e/MCL_project/chipseq_h3_ubf/"
path_to_raw_fastq="$path_to_project""raw_data/"
path_to_clean_fastq="$path_to_project""raw_data/clean_fq/"
path_to_ref="/mnt/e/MCL_project/chipseq_nucleolin/ref"
path_to_bwa_files="$path_to_project""bwa/"
path_to_picard="/mnt/e/MCL_project/"
path_to_bbmap="/mnt/e/MCL_project/bbmap/"
#path_to_macs2_files="$path_to_project""macs2_pe/"
#path_to_macs2_files_withdup="$path_to_project""macs2_pe_dup/"
#path_to_macs2_files_lessstringent="$path_to_project""macs2_lessstringent/"
#path_to_idr="$path_to_project""idr/"
#path_to_bigwig="$path_to_project""bigwigs/"
genome_fa="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
samples=(L1a L2a L3a L1b L2b L3b L1c L2c L3c Linput G1a G2a G3a G1b G2b G3b G1c G2c G3c Ginput)

BEGINCOMMENT 
### FastQC quality analysis
for i in ${samples[*]};do
    fastqc -o "$path_to_raw_fastq""fastqc/" -t 12 "$path_to_raw_fastq""${i}"/*.fq.gz &
done
ENDCOMMENT

### Adapter trimming 
for i in ${samples[*]};do
    "$path_to_bbmap"bbduk.sh in1="$path_to_raw_fastq""${i}"/${i}_1.fq.gz in2="$path_to_raw_fastq""${i}"/${i}_2.fq.gz \
    out1="$path_to_clean_fastq"${i}_clean_1.fq out2="$path_to_clean_fastq"${i}_clean_2.fq \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo &
done
