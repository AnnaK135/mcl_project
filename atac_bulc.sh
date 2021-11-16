#!/bin/bash


### For commenting out pieces of the code ###
[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

path_to_project="/mnt/e/MCL_project/atac/"
path_to_tools="/mnt/e/tools/"
path_to_ref="/mnt/e/MCL_project/chipseq_nucleolin/ref"
path_to_raw_fastq="$path_to_project""fastq/"
path_to_clean_fastq="$path_to_project""fastq_clean/"
path_to_fastqc="$path_to_project""fastqc/"
path_to_bwa_files="$path_to_project""bwa/"
genome_fa="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
index="GRCh38_noalt_as.1.bt2"
srr=(SRR891268 SRR891269 SRR891270 SRR891271)

BEGINCOMMENT

### FastQC quality analysis
fastqc "$path_to_raw_fastq"* -o "$path_to_fastqc" &

### Removing adapters with NGmerge
for i in ${srr[*]};do
    "$path_to_tools"NGmerge/NGmerge -a -1 "$path_to_raw_fastq"${i}_1.fastq -2 "$path_to_raw_fastq"${i}_2.fastq -o "$path_to_clean_fastq"${i} -v -n 10 -u 41 &
done
ENDCOMMENT

### FastQC quality analysis of clean reads
fastqc "$path_to_clean_fastq"* -o "$path_to_fastqc" &

### Alignment with Bowtie, stats and sorting
for i in ${srr[*]};do
    bowtie2 --very-sensitive -x "$path_to_ref""$index" -k 10 -1 "$path_to_clean_fastq"${i}_1.fastq - 2 "$path_to_clean_fastq"${i}_clean_2.fastq -p 10 \
        | samtools view -u - \
        | samtools sort -n - o "$path_to_bwa_files"${i}_sorted.bam -
    samtools flagstat "$path_to_bwa_files"${i}_sorted.bam > "$path_to_bwa_files"${i}.mapping_stat
done

BEGINCOMMENT
### Peak calling with Genrich
#(in the ATACseq mode -j, removing duplicates -r and mitochondrial reads -e chrM and for multiple replicates at once)
"$path_to_tools"Genrich/Genrich -j -r -e chrM

ENDCOMMENT