#!/bin/bash

### For commenting out pieces of the code ###
[ -z $BASH ] || shopt -s expand_aliases
alias BEGINCOMMENT="if [ ]; then"
alias ENDCOMMENT="fi"

### Introducing variables for paths so we can change them later and reuse the script ###
path_to_project="/mnt/e/MCL_project/chipseq_h3_ubf/"
path_to_raw_fastq="$path_to_project""raw_data/"
path_to_clean_fastq="$path_to_project""raw_data/clean_fq/" 
path_to_ref="/mnt/e/MCL_project/ref/"
path_to_bwa_files="$path_to_project""bwa/"
path_to_picard="/mnt/e/MCL_project/"
path_to_bbmap="/mnt/e/MCL_project/bbmap/"
path_to_macs2_files="$path_to_project""macs2/"
path_to_idr="$path_to_project""idr/"
path_to_chipr="$path_to_project""chipr/"
path_to_bigwig="$path_to_project""bigwigs/"
genome_fa="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
rose.lcl=(L1a L2a L3a L1b L2b L3b L1c L2c L3c Linput)
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

for i in ${broads_l[*]};do
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

for i in ${broads_l[*]};do
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

for i in ${broads_l[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_filtered.bam \
                   -c "$path_to_bwa_files"Linput_sorted_marked_filtered.bam \  
                   -n ${i}_lessstringent_broad -f BAMPE -g hs --pvalue 1e-3 --fix-bimodal --broad \
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_lesstringent_broad_macs2.log 
done

ENDCOMMENT

### Calling broad peaks (H3K27Ac)

broads_la=(L1b L2b L3b Linput)
broads_ga=(G1b G2b G3b Ginput)

#no duplicates
for i in ${broads_ga[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_filtered.bam \
                   -c "$path_to_bwa_files"Ginput_sorted_marked_filtered.bam \  
                   -n ${i}_nodup_broad -f BAMPE -g hs --fix-bimodal --broad \
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_nodup_broad_macs2.log 
done

for i in ${broads_la[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_duplicates.bam \
                   -c "$path_to_bwa_files"Linput_sorted_marked_duplicates.bam \  
                   -n ${i}_nodup_broad -f BAMPE -g hs --fix-bimodal --broad \
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_nodup_broad_macs2.log 
done

#duplicates
for i in ${broads_ga[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_duplicates.bam \
                   -c "$path_to_bwa_files"Ginput_sorted_marked_duplicates.bam \  
                   -n ${i}_dup_broad -f BAMPE -g hs --fix-bimodal --broad \
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_dup_broad_macs2.log 
done

for i in ${broads_la[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_duplicates.bam \
                   -c "$path_to_bwa_files"Linput_sorted_marked_duplicates.bam \  
                   -n ${i}_dup_broad -f BAMPE -g hs --fix-bimodel --broad \
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_dup_broad_macs2.log
done

#less stringent
for i in ${broads_ga[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_filtered.bam \
                   -c "$path_to_bwa_files"Ginput_sorted_marked_filtered.bam \  
                   -n ${i}_lessstringent_broad -f BAMPE -g hs --pvalue 1e-3 --fix-bimodal --broad \
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_lesstringent_broad_macs2.log 
done

for i in ${broads_la[*]};do
    macs2 callpeak -t "$path_to_bwa_files"${i}_sorted_marked_filtered.bam \
                   -c "$path_to_bwa_files"Linput_sorted_marked_filtered.bam \  
                   -n ${i}_lessstringent_broad -f BAMPE -g hs --pvalue 1e-3 --fix-bimodal --broad \
                   --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_lesstringent_broad_macs2.log 
done

BEGINCOMMENT

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

### Merging replicates with IDR ### 
## for narrow peaks (UBF, H3K27Ac)

idr --samples "$path_to_macs2_files"L1a_lessstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"L2a_lessstringent_peaks_sorted.narrowPeak \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file "$path_to_idr"LCL_UBF12_idr \
    --plot \
    --log-output-file "$path_to_idr"LCL_UBF12_idr.log

idr --samples "$path_to_macs2_files"G1a_lessstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"G2a_lessstringent_peaks_sorted.narrowPeak \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file "$path_to_idr"GRANTA_UBF12_idr \
    --plot \
    --log-output-file "$path_to_idr"GRANTA_UBF12_idr.log

idr --samples "$path_to_macs2_files"L1b_lessstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"L2b_lessstringent_peaks_sorted.narrowPeak \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file "$path_to_idr"LCL_H3K27Ac12_idr \
    --plot \
    --log-output-file "$path_to_idr"LCL_H3K27Ac12_idr.log

idr --samples "$path_to_macs2_files"G1b_lessstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"G2b_lessstringent_peaks_sorted.narrowPeak \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file "$path_to_idr"GRANTA_H3K27Ac12_idr \
    --plot \
    --log-output-file "$path_to_idr"GRANTA_H3K27Ac12_idr.log

### Sorting non-stringent broadPeaks by p-value ###
for i in ${samples[*]};do    
    sort -k8,8nr "$path_to_macs2_files"${i}_lessstringent_broad_peaks.broadPeak > "$path_to_macs2_files"${i}_lessstringent_peaks_sorted.broadPeak &
done

### Merging broad peaks with IDR
idr --samples "$path_to_macs2_files"L1c_lessstringent_peaks_sorted.broadPeak "$path_to_macs2_files"L2c_lessstringent_peaks_sorted.broadPeak \
    --input-file-type broadPeak \
    --rank p.value \
    --output-file "$path_to_idr"LCL_H3K36me12_idr \
    --plot \
    --log-output-file "$path_to_idr"LCL_H3K36me12_idr.log

idr --samples "$path_to_macs2_files"G1c_lessstringent_peaks_sorted.broadPeak "$path_to_macs2_files"G3c_lessstringent_peaks_sorted.broadPeak \
    --input-file-type broadPeak \
    --rank p.value \
    --output-file "$path_to_idr"GRANTA_H3K36me13_idr \
    --plot \
    --log-output-file "$path_to_idr"GRANTA_H3K36me13_idr.log

### Merging replicates with ChIP-R

chipr -i "$path_to_macs2_files"L1a_lesstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"L2a_lesstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"L3a_lesstringent_peaks_sorted.narrowPeak \
        -m 2 -o "$_path_to_chipr"LCL_UBF_chipr &

chipr -i "$path_to_macs2_files"G1a_lesstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"G2a_lesstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"G3a_lesstringent_peaks_sorted.narrowPeak \
        -m 2 -o "$_path_to_chipr"GRANTA_UBF_chipr &

chipr -i "$path_to_macs2_files"L1b_lesstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"L2b_lesstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"L3b_lesstringent_peaks_sorted.narrowPeak \
        -m 2 -o "$_path_to_chipr"LCL_H3K27Ac_chipr &

chipr -i "$path_to_macs2_files"G1b_lesstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"G2b_lesstringent_peaks_sorted.narrowPeak "$path_to_macs2_files"G3b_lesstringent_peaks_sorted.narrowPeak \
        -m 2 -o "$_path_to_chipr"GRANTA_H3K27Ac_chipr &

chipr -i "$path_to_macs2_files"L1c_lessstringent_peaks_sorted.broadPeak "$path_to_macs2_files"L2c_lessstringent_peaks_sorted.broadPeak "$path_to_macs2_files"L3c_lessstringent_peaks_sorted.broadPeak \
        -m 2 -o "$_path_to_chipr"LCL_H3K36me13_chipr &

chipr -i "$path_to_macs2_files"G1c_lessstringent_peaks.narrowPeak "$path_to_macs2_files"G2c_lessstringent_peaks.narrowPeak "$path_to_macs2_files"G3c_lessstringent_peaks.narrowPeak \
        -m 2 -o "$_path_to_chipr"GRANTA_H3K36me13_chipr &

### Creating score matrices for chr11, chr14 and all chromosomes. UBF
computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R "$path_to_ref"refseq_chr11_genes.bed \
    -S "$path_to_bigwig"ubf/*.bw \
    --skipZeros -p 6 \
    -o "$path_to_project"matrix/matrix_ubf_chr11_tss.gz
    --outFileSortedRegions "$path_to_project"matrix/ubf_chr11_tss.bed &

computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R "$path_to_ref"refseq_chr14_genes.bed \
    -S "$path_to_bigwig"ubf/*.bw \
    --skipZeros -p 6 \
    -o "$path_to_project"matrix/matrix_ubf_chr14_tss.gz
    --outFileSortedRegions "$path_to_project"matrix/ubf_chr14_tss.bed &

computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R "$path_to_ref"refseq_all_genes.bed \
    -S "$path_to_bigwig"ubf/*.bw \
    --skipZeros -p 6 \
    -o "$path_to_project"matrix/matrix_ubf_all_tss.gz
    --outFileSortedRegions "$path_to_project"matrix/ubf_all_tss.bed &


### Creating score matrices for chr11, chr14 and all chromosomes. H3K27Ac

computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R "$path_to_ref"refseq_chr11_genes.bed \
    -S "$path_to_bigwig"h3k27ac/*.bw \
    --skipZeros -p 6 \
    -o "$path_to_project"matrix/matrix_h3k27ac_chr11_tss.gz
    --outFileSortedRegions "$path_to_project"matrix/h3k27ac_chr11_tss.bed &

computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R "$path_to_ref"refseq_chr14_genes.bed \
    -S "$path_to_bigwig"h3k27ac/*.bw \
    --skipZeros -p 6 \
    -o "$path_to_project"matrix/matrix_h3k27ac_chr14_tss.gz
    --outFileSortedRegions "$path_to_project"matrix/h3k27ac_chr14_tss.bed &

computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R "$path_to_ref"refseq_all_genes.bed \
    -S "$path_to_bigwig"h3k27ac/*.bw \
    --skipZeros -p 6 \
    -o "$path_to_project"matrix/matrix_h3k27ac_all_tss.gz
    --outFileSortedRegions "$path_to_project"matrix/h3k27ac_all_tss.bed &

### Creating score matrices for chr11, chr14 and all chromosomes. H3K36me3

computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R "$path_to_ref"refseq_chr11_genes.bed \
    -S "$path_to_bigwig"h3k36me3/*.bw \
    --skipZeros -p 6 \
    -o "$path_to_project"matrix/matrix_h3k36me3_chr11_tss.gz
    --outFileSortedRegions "$path_to_project"matrix/h3k36me3_chr11_tss.bed &

computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R "$path_to_ref"refseq_chr14_genes.bed \
    -S "$path_to_bigwig"h3k36me3/*.bw \
    --skipZeros -p 6 \
    -o "$path_to_project"matrix/matrix_h3k36me3_chr14_tss.gz
    --outFileSortedRegions "$path_to_project"matrix/h3k36me3_chr14_tss.bed &

computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R "$path_to_ref"refseq_all_genes.bed \
    -S "$path_to_bigwig"h3k36me3/*.bw \
    --skipZeros -p 6 \
    -o "$path_to_project"matrix/matrix_h3k36me3_all_tss.gz
    --outFileSortedRegions "$path_to_project"matrix/h3k36me3_all_tss.bed &


ENDCOMMENT

