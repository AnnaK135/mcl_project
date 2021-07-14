#!/bin/bash

### For commenting out pieces of the code ###
[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

### Introducing variables for paths so we can change them later and reuse the script ###
path_to_project="/mnt/e/MCL_project/chipseq_nucleolin/"
path_to_raw_fastq="$path_to_project""raw_data/CleanFq/"
path_to_ref="$path_to_project""ref/"
path_to_bwa_files="$path_to_project""bwa/"
path_to_picard="/mnt/e/MCL_project/"
path_to_macs2_files="$path_to_project""macs2/"
genome_fa="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

BEGINCOMMENT
#### Dependencies ###
#sudo apt-get update
#sudo apt install bwa #version 0.7.17 at the moment of execution
#sudo apt install samtools
#pip3 install macs2

#### Getting the fastq prefixes ### 
ls "$path_to_raw_fastq" | sed -e 's/\_[1,2].fq.gz$//' | uniq > "$path_to_project"fastq_prefixes

#### BWA mapping and sorting BAMs with Samtools ###
#(BWA index was pre-downloaded from usca and stored in ref)

for i in $(cat "$path_to_project"fastq_prefixes);do
    bwa mem -t 20 -M "$path_to_ref""$genome_fa" "$path_to_raw_fastq"${i}_1.fq.gz "$path_to_raw_fastq"${i}_2.fq.gz | samtools view -bS > "$path_to_bwa_files"${i}.bam
    echo "Alignment for "${i}" completed! BAM file is saved to" "$path_to_bwa_files"
    samtools flagstat "$path_to_bwa_files"${i}.bam > "$path_to_bwa_files"${i}.mapping_stat
    echo "Stats for "${i}" saved to" ""$path_to_bwa_files"${i}.mapping_stat"
    echo "Sorting "${i}".bam ..." 
    samtools sort "$path_to_bwa_files"${i}.bam > "$path_to_bwa_files"${i}_sorted.bam
    echo "Sorting done"      
done
# -M mark shorter split hits as secondary (for Picard compatibility)

#### Marking duplicates with Picard ###
for i in $(cat "$path_to_project"fastq_prefixes);do
    java -jar "$path_to_picard"/picard.jar MarkDuplicates \
    I="$path_to_bwa_files"${i}_sorted.bam \
    O="$path_to_bwa_files"${i}_sorted_marked_duplicates.bam\
    M="$path_to_bwa_files"${i}_sorted_dup_metrics.txt \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT
    echo "Marking for "${i}" completed! Marked BAM file is saved to" "$path_to_bwa_files"
done

### Filtering (and indexing for the genome browser) ###

#filtering unmapped reads and reads (flag 4) with mapping quality < 5 and Picard-Marked duplicates (flag 1024)
for i in $(cat "$path_to_project"fastq_prefixes);do
    samtools index -b -@ 20 "$path_to_bwa_files"${i}_sorted_marked_duplicates.bam 
    echo "Indexing for "${i}" completed! BAI file is saved to" "$path_to_bwa_files"
    samtools view -b -F 4 -q 5 -@ 20 "$path_to_bwa_files"${i}_sorted_marked_duplicates.bam | samtools view -b -F 1024 -@ 20 > "$path_to_bwa_files"${i}_sorted_marked_filtered.bam
    echo "Filtering for "${i}" completed! Filtered BAM is saved to" "$path_to_bwa_files"
done

for i in $(cat "$path_to_project"fastq_prefixes);do
    run_spp.R -c="$path_to_bwa_files"${i}_sorted_marked_filtered.bam -rf -out="$path_to_bwa_files"${i}_sorted_marked_filtered.phantom.out
    echo "Fragment size and quality parameters estimated for "${i}" ! Output saved to" ""$path_to_bwa_files"${i}_sorted_marked_filtered.phantom.out"
done

ENDCOMMENT

path_to_macs2_files="$path_to_project""macs2/"
for i in $(cat "$path_to_project"short_prefixes);do
    macs2 callpeak -t "$path_to_bwa_files"${i}IP_clean_sorted_marked_filtered.bam \
                    -c "$path_to_bwa_files"${i}IN_clean_sorted_marked_filtered.bam \
                    -n ${i} -f BAM -g hs \
                    --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_macs2.log &
done
