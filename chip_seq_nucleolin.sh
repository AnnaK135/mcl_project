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
path_to_macs2_files="$path_to_project""macs2_pe/"
path_to_macs2_files_withdup="$path_to_project""macs2_pe_dup/"
path_to_macs2_files_lessstringent="$path_to_project""macs2_lessstringent/"
path_to_idr="$path_to_project""idr/"
path_to_bigwig="$path_to_project""bigwigs_bmp/"
genome_fa="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
samples=(LCL GRANTA)

BEGINCOMMENT
#### Dependencies ###
#sudo apt-get update
#sudo apt install bwa #version 0.7.17 at the moment of execution
#sudo apt install samtools
#pip3 install macs2
#https://github.com/nboley/idr.git - IDR package (to be installed independently from GitHub, version 2.0.3)

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

### MACS2 peak calling with removed duplicates ###

path_to_macs2_files="$path_to_project""macs2_pe/"
for i in $(cat "$path_to_project"short_prefixes);do
    macs2 callpeak -t "$path_to_bwa_files"${i}IP_clean_sorted_marked_filtered.bam \
                    -c "$path_to_bwa_files"${i}IN_clean_sorted_marked_filtered.bam \
                    -n ${i} -g hs -f BAMPE \
                    --outdir "$path_to_macs2_files" 2> "$path_to_macs2_files"${i}_macs2.log &
done

### MACS2 Peak calling with duplicates (for differerntial binding analysis) ###

for i in $(cat "$path_to_project"short_prefixes);do
    macs2 callpeak -t "$path_to_bwa_files"${i}IP_clean_sorted.bam \
                   -c "$path_to_bwa_files"${i}IN_clean_sorted.bam \
                   -n ${i} -g hs -f BAMPE --keep-dup all\
                   --outdir "$path_to_macs2_files_withdup" 2> "$path_to_macs2_files_withdup"${i}_macs2.log &
done

### Less strigent MACS2 peak calling for merging replicates with IDR ###

for i in $(cat "$path_to_project"short_prefixes);do
    macs2 callpeak -t "$path_to_bwa_files"${i}IP_clean_sorted_marked_filtered.bam \
                   -c "$path_to_bwa_files"${i}IN_clean_sorted_marked_filtered.bam \  
                   -n ${i} -f BAMPE -g hs --pvalue 1e-3 \
                   --outdir "$path_to_macs2_files_lessstringent" 2> "$path_to_macs2_files_lessstringent"${i}_macs2.log &
done

### Sorting non-stringent narrowPeaks by p-value ###

for i in $(cat "$path_to_project"short_prefixes);do    
    sort -k8,8nr "$path_to_macs2_files_lessstringent"${i}_peaks.narrowPeak > "$path_to_macs2_files_lessstringent"${i}_peaks_sorted.narrowPeak &
done

### Merging replicates with IDR ### 

for i in ${samples[*]};do
    idr --samples "$path_to_macs2_files_lessstringent"${i}a_peaks_sorted.narrowPeak "$path_to_macs2_files_lessstringent"${i}b_peaks_sorted.narrowPeak \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file "$path_to_idr"${i}_idr \
    --plot \
    --log-output-file "$path_to_idr"${i}_idr.log &
done

### Creating pseudoreplicates for IP ### ---RETRY FROM HERE TILL THE MARK
for i in ${samples[*]};do
    samtools merge -u "$path_to_bwa_files"${i}_mergedIP.bam "$path_to_bwa_files"${i}aIP_clean_sorted_marked_filtered.bam "$path_to_bwa_files"${i}bIP_clean_sorted_marked_filtered.bam
    samtools view -H "$path_to_bwa_files"${i}_mergedIP.bam > "$path_to_bwa_files"${i}_mergedIP_header.sam &&
    nlines=$(samtools view "$path_to_bwa_files"${i}_merged.bam | wc -l ) && # Number of reads in the merged BAM 
    nlines=$(( (nlines + 1) / 2 )) &&  # half that number 
    samtools view "$path_to_bwa_files"${i}_mergedIP.bam | shuf - | split -d -l ${nlines} "$path_to_bwa_files"${i}_pseudoIP &&
    cat "$path_to_bwa_files"${i}_mergedIP_header.sam "$path_to_bwa_files"${i}_pseudoIP00 | samtools view -bS - > "$path_to_bwa_files"${i}_pseudoIP00.bam
    cat "$path_to_bwa_files"${i}_mergedIP_header.sam "$path_to_bwa_files"${i}_pseudoIP01 | samtools view -bS - > "$path_to_bwa_files"${i}_pseudoIP01.bam
done

### Creating pseudoreplicates for IN ###
for i in ${samples[*]};do
    samtools merge -u "$path_to_bwa_files"${i}_mergedIN.bam "$path_to_bwa_files"${i}aIN_clean_sorted_marked_filtered.bam "$path_to_bwa_files"${i}bIN_clean_sorted_marked_filtered.bam
    samtools view -H "$path_to_bwa_files"${i}_mergedIN.bam > "$path_to_bwa_files"${i}_mergedIN_header.sam &&
    nlines=$(samtools view "$path_to_bwa_files"${i}_merged.bam | wc -l ) && # Number of reads in the merged BAM 
    nlines=$(( (nlines + 1) / 2 )) &&  # half that number 
    samtools view "$path_to_bwa_files"${i}_mergedIN.bam | shuf - | split -d -l ${nlines} "$path_to_bwa_files"${i}_pseudoIN &&
    cat "$path_to_bwa_files"${i}_mergedIN_header.sam "$path_to_bwa_files"${i}_pseudoIN00 | samtools view -bS - > "$path_to_bwa_files"${i}_pseudoIN00.bam
    cat "$path_to_bwa_files"${i}_mergedIN_header.sam "$path_to_bwa_files"${i}_pseudoIN01 | samtools view -bS - > "$path_to_bwa_files"${i}_pseudoIN01.bam
done

### Peak calling and sorting on pseudoreplicates ###

for i in ${samples[*]};do
    for k in 00 01;do
    macs2 callpeak -t "$path_to_bwa_files"${i}_pseudoIP"$k".bam \
                   -c "$path_to_bwa_files"${i}_pseudoIN"$k".bam \  
                   -n ${i}_pseudo"$k" -f BAMPE -g hs -p 1e-3 \
                   --outdir "$path_to_macs2_files_lessstringent" 2> "$path_to_macs2_files_lessstringent"${i}_pseudo"$k"_macs2.log &
    done
done

for i in ${samples[*]};do
    for k in 00 01;do
    sort -k8,8nr "$path_to_macs2_files_lessstringent"${i}_pseudo"$k"_peaks.narrowPeak > "$path_to_macs2_files_lessstringent"${i}_pseudo"$k"_peaks_sorted.narrowPeak &
    done
done

### IDR on pseudoreplicates ###

for i in ${samples[*]};do
idr --samples "$path_to_macs2_files_lessstringent"${i}_pseudo00_peaks_sorted.narrowPeak "$path_to_macs2_files_lessstringent"${i}${i}_pseudo01_peaks_sorted.narrowPeak \
    --input-file-type narrowPeak \
    --rank p.value \
    --output-file "$path_to_idr"${i}_pseudo_idr \
    --plot \
    --log-output-file "$path_to_idr"${i}_pseudo_idr.log &
done

###---THE MARK

### Creating BigWig files for visualistaion ### 

### bamCompare for each IP-IN pair (normalisation of IP to input)
for i in $(cat "$path_to_project"short_prefixes);do
    bamCompare -b1 "$path_to_bwa_files"${i}IP_clean_sorted_marked_duplicates.bam \
               -b2 "$path_to_bwa_files"${i}IN_clean_sorted_marked_duplicates.bam \
               -o "$path_to_bigwig"${i}_log2ratio.bw \
               --binSize 10 --normalizeUsing BPM --smoothLength 30 \
	       --scaleFactorsMethod None \
               --extendReads 150 --centerReads -p 4 -v 2> "$path_to_bigwig"${i}_log2ratio.log &
done

### bamCoverage for each file (normalization to sequencing depth)
for i in $(cat "$path_to_project"short_prefixes);do
    bamCoverage -b "$path_to_bwa_files"${i}IP_clean_sorted_marked_duplicates.bam \
                -o "$path_to_bigwig"${i}IP_clean_sorted_marked_duplicates_coverage.bw \
                --binSize 10 --normalizeUsing BPM --smoothLength 30 \
                --extendReads 150 --centerReads -p 4 -v 2> "$path_to_bigwig"${i}IP_coverage.log &
done

for i in $(cat "$path_to_project"short_prefixes);do
    bamCoverage -b "$path_to_bwa_files"${i}IN_clean_sorted_marked_duplicates.bam \
                -o "$path_to_bigwig"${i}IN_clean_sorted_marked_duplicates_coverage.bw \
                --binSize 10 --normalizeUsing BPM --smoothLength 30 \
                --extendReads 150 --centerReads -p 4 -v 2> "$path_to_bigwig"${i}IN_coverage.log &
done

ENDCOMMENT

### Creating score matrices for chr11, chr14 and all chromosomes.

computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R "$path_to_ref"refseq_chr11_genes.bed \
    -S "$path_to_bigwig"*.bw \
    --skipZeros -p 6 \
    -o "$path_to_project"matrix/matrix_nuc_chr11_tss.gz
    --outFileSortedRegions "$path_to_project"matrix/nuc_chr11_tss.bed

computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R "$path_to_ref"refseq_chr14_genes.bed \
    -S "$path_to_bigwig"*.bw \
    --skipZeros -p 6 \
    -o "$path_to_project"matrix/matrix_nuc_chr14_tss.gz
    --outFileSortedRegions "$path_to_project"matrix/nuc_chr14_tss.bed 

computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R "$path_to_ref"refseq_all_genes.bed \
    -S "$path_to_bigwig"*.bw \
    --skipZeros -p 6 \
    -o "$path_to_project"matrix/matrix_nuc_all_tss.gz
    --outFileSortedRegions "$path_to_project"matrix/nuc_all_tss.bed 
