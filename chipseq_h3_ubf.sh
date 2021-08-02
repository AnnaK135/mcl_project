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
path_to_macs2_files="$path_to_project""macs2/"
path_to_idr="$path_to_project""idr/"
path_to_bigwig="$path_to_project""bigwigs/"
genome_fa="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
samples=(L1a L2a L3a L1b L2b L3b L1c L2c L3c Linput G1a G2a G3a G1b G2b G3b G1c G2c G3c Ginput)
lcl=(L1a L2a L3a L1b L2b L3b L1c L2c L3c Linput)
granta=(G1a G2a G3a G1b G2b G3b G1c G2c G3c Ginput)

BEGINCOMMENT
### FastQC quality analysis
for i in ${samples[*]};do
    fastqc -o "$path_to_raw_fastq""fastqc/" -t 12 "$path_to_raw_fastq""${i}"/*.fq.gz &
done
 
### Adapter trimming

for i in ${samples[*]};do
    "$path_to_bbmap"bbduk.sh in1="$path_to_raw_fastq""${i}"/${i}_1.fq.gz in2="$path_to_raw_fastq""${i}"/${i}_2.fq.gz \
    out1="$path_to_clean_fastq"${i}_clean_1.fq out2="$path_to_clean_fastq"${i}_clean_2.fq \
    ref="$path_to_bbmap"resources/adapters.fa \
    ktrim=r k=23 mink=11 hdist=1 tpe tbo
done

### Alignment with BWA-MEM, stats and sorting
for i in ${samples[*]};do
    bwa mem -t 20 -M "$path_to_ref"$genome_fa "$path_to_clean_fastq"${i}_clean_1.fq "$path_to_clean_fastq"${i}_clean_2.fq | samtools view -bS > "$path_to_bwa_files"${i}.bam
    samtools flagstat "$path_to_bwa_files"${i}.bam > "$path_to_bwa_files"${i}.mapping_stat
    samtools sort "$path_to_bwa_files"${i}.bam > "$path_to_bwa_files"${i}_sorted.bam
done

#### Marking duplicates with Picard ###
for i in ${samples[*]};do
    java -jar "$path_to_picard"/picard.jar MarkDuplicates \
    I="$path_to_bwa_files"${i}_sorted.bam \
    O="$path_to_bwa_files"${i}_sorted_marked_duplicates.bam\
    M="$path_to_bwa_files"${i}_sorted_dup_metrics.txt \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT
done

### Filtering (and indexing for the genome browser) ###
#filtering unmapped reads and reads (flag 4) with mapping quality < 5 and Picard-Marked duplicates (flag 1024)
for i in ${samples[*]};do
    samtools index -b -@ 20 "$path_to_bwa_files"${i}_sorted_marked_duplicates.bam 
    samtools view -b -F 4 -q 5 -@ 20 "$path_to_bwa_files"${i}_sorted_marked_duplicates.bam | samtools view -b -F 1024 -@ 20 > "$path_to_bwa_files"${i}_sorted_marked_filtered.bam
done


### Calling narrow peaks (UBF and H3K27Ac, and H3K36me3 to try)

### MACS2 peak calling with removed duplicates ###
for i in ${lcl[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_filtered.bam \
                   -c "$path_to_bwa_files"Linput_sorted_marked_filtered.bam \
                   --fix-bimodal -f BAMPE -n ${i}_nodup -g hs \
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_nodup_macs2.log 
done

for i in ${granta[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_filtered.bam \
                   -c "$path_to_bwa_files"Ginput_sorted_marked_filtered.bam \
                   -n ${i}_nodup -g hs -f BAMPE --fix-bimodal\
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_nodup_macs2.log 
done

### MACS2 Peak calling with duplicates (for differerntial binding analysis) ###
for i in ${lcl[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_duplicates.bam \
                   -c "$path_to_bwa_files"Linput_sorted_marked_duplicates.bam \
                   -n ${i}_dup -g hs -f BAMPE --keep-dup all --fix-bimodal\
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_dup_macs2.log 
done

for i in ${granta[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_duplicates.bam \
                   -c "$path_to_bwa_files"Ginput_sorted_marked_duplicates.bam \
                   -n ${i}_dup -g hs -f BAMPE --keep-dup all --fix-bimodal\
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_dup_macs2.log 
done

### Less strigent MACS2 peak calling for merging replicates with IDR ###
for i in ${lcl[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_filtered.bam \
                   -c "$path_to_bwa_files"Linput_sorted_marked_filtered.bam \
                   -n ${i}_lessstringent -f BAMPE -g hs --pvalue 1e-3 --fix-bimodal \ 
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_lesstringent_macs2.log &
done

for i in ${granta[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_filtered.bam \
                   -c "$path_to_bwa_files"Ginput_sorted_marked_filtered.bam \  
                   -n ${i}_lessstringent -f BAMPE -g hs --pvalue 1e-3 --fix-bimodal\
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_lesstringent_macs2.log 
done

### Calling broad peaks (H3K36me3)

broads_l=(L1c L2c L3c Linput)
broads_g=(G1c G2c G3c Ginput)

#no duplicates
for i in ${broads_g[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_filtered.bam \
                   -c "$path_to_bwa_files"Ginput_sorted_marked_filtered.bam \  
                   -n ${i}_nodup_broad -f BAMPE -g hs --fix-bimodal --broad \
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_nodup_broad_macs2.log 
done

for i in ${lcl_g[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_duplicates.bam \
                   -c "$path_to_bwa_files"Linput_sorted_marked_duplicates.bam \  
                   -n ${i}_nodup_broad -f BAMPE -g hs --fix-bimodal --broad \
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_nodup_broad_macs2.log 
done

#duplicates
for i in ${broads_g[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_duplicates.bam \
                   -c "$path_to_bwa_files"Ginput_sorted_marked_duplicates.bam \  
                   -n ${i}_dup_broad -f BAMPE -g hs --fix-bimodal --broad \
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_dup_broad_macs2.log 
done

for i in ${lcl_g[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_duplicates.bam \
                   -c "$path_to_bwa_files"Linput_sorted_marked_duplicates.bam \  
                   -n ${i}_dup_broad -f BAMPE -g hs --fix-bimodel --broad \
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_dup_broad_macs2.log
done

#less stringent
for i in ${broads_g[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_filtered.bam \
                   -c "$path_to_bwa_files"Ginput_sorted_marked_filtered.bam \  
                   -n ${i}_lessstringent_broad -f BAMPE -g hs --pvalue 1e-3 --fix-bimodal --broad \
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_lesstringent_broad_macs2.log 
done

for i in ${lcl_g[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_filtered.bam \
                   -c "$path_to_bwa_files"Linput_sorted_marked_filtered.bam \  
                   -n ${i}_lessstringent_broad -f BAMPE -g hs --pvalue 1e-3 --fix-bimodal --broad \
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_lesstringent_broad_macs2.log 
done

### Sorting non-stringent narrowPeaks by p-value ###

for i in ${samples[*]};do
    sort -k8,8nr "$path_to_macs2_files"${i}_lessstringent_peaks.narrowPeak > "$path_to_macs2_files"${i}_lesstringent_peaks_sorted.narrowPeak
done

### bamCompare for each IP-IN pair (normalisation of IP to input)
for i in ${lcl[*]};do
    bamCompare -b1 "$path_to_bwa_files"${i}_sorted_marked_duplicates.bam \
               -b2 "$path_to_bwa_files"Linput_sorted_marked_duplicates.bam \
               -o "$path_to_bigwig"${i}_log2ratio.bw \
               --binSize 10 --normalizeUsing BPM --smoothLength 30 \
               --scaleFactorsMethod None \
               --extendReads 150 --centerReads -p 4 -v 2> "$path_to_bigwig"${i}_log2ratio.log &
done

for i in ${granta[*]};do
    bamCompare -b1 "$path_to_bwa_files"${i}_sorted_marked_duplicates.bam \
               -b2 "$path_to_bwa_files"Ginput_sorted_marked_duplicates.bam \
               -o "$path_to_bigwig"${i}_log2ratio.bw \
               --binSize 10 --normalizeUsing BPM --smoothLength 30 \
               --scaleFactorsMethod None \
               --extendReads 150 --centerReads -p 4 -v 2> "$path_to_bigwig"${i}_log2ratio.log &
done

### bamCoverage for each file (normalization to sequencing depth)
for i in ${samples[*]};do
    bamCoverage -b "$path_to_bwa_files"${i}_sorted_marked_duplicates.bam \
                -o "$path_to_bigwig"${i}_sorted_marked_duplicates_coverage.bw \
                --binSize 10 --normalizeUsing BPM --smoothLength 30 \
                --extendReads 150 --centerReads -p 4 -v 2> "$path_to_bigwig"${i}_coverage.log &
done

ENDCOMMENT

### Merging replicates with IDR ### 
## for narrow peaks (UBF, H3K27Ac)

idr --samples "$path_to_macs2_files"L1a_lessstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"L2a_lessstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"L3a_lessstringent_peaks_sorted.narrowPeak \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file "$path_to_idr"LCL_UBF_idr \
    --plot \
    --log-output-file "$path_to_idr"LCL_UBF_idr.log

idr --samples "$path_to_macs2_files"G1a_lessstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"G2a_lessstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"G3a_lessstringent_peaks_sorted.narrowPeak \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file "$path_to_idr"GRANTA_UBF_idr \
    --plot \
    --log-output-file "$path_to_idr"GRANTA_UBF_idr.log

idr --samples "$path_to_macs2_files"L1b_lessstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"L2b_lessstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"L3b_lessstringent_peaks_sorted.narrowPeak \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file "$path_to_idr"LCL_H3K27Ac_idr \
    --plot \
    --log-output-file "$path_to_idr"LCL_H3K27Ac_idr.log

idr --samples "$path_to_macs2_files"G1b_lessstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"G2b_lessstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"G3b_lessstringent_peaks_sorted.narrowPeak \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file "$path_to_idr"GRANTA_H3K27Ac_idr \
    --plot \
    --log-output-file "$path_to_idr"GRANTA_H3K27Ac_idr.log
