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
path_to_diffbind="$path_to_project""diffbind_h3k27ac/results/"
path_to_bigwig="$path_to_project""bigwigs/h3k27ac/"
path_to_deseq="$path_to_project""deseq/"

samples=(L1b L2b L3b Linput G1b G2b G3b Ginput)
lcl=(L1b L2b L3b)
granta=(G1b G2b G3b)

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

chipr -i "$path_to_rose_files"rose_output_G1b/GRANTA_H3K27Ac_chipr_all_Gateway_SuperEnhancers.bed \
         "$path_to_rose_files"rose_output_G2b/GRANTA_H3K27Ac_chipr_all_Gateway_SuperEnhancers.bed \
         "$path_to_rose_files"rose_output_G3b/GRANTA_H3K27Ac_chipr_all_Gateway_SuperEnhancers.bed \
          -m 2 -o "$path_to_rose_files"GRANTA_SEs_chipr 


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

### Computing the matrix for the differentially h3k27ac enriched regions  
computeMatrix reference-point --referencePoint center \
   -b 3000 -a 3000 \
   -R "$path_to_diffbind"diffenriched_h3k27ac.bed \
   -S "$path_to_bigwig"*.bw \
   --skipZeros -p 6 \
   -o "$path_to_project"matrix/matrix_diffenriched_h3k27ac_centered.gz \
   --outFileSortedRegions "$path_to_project"matrix/matrix_diffenriched_h3k27ac_centered.bed
 
computeMatrix reference-point --referencePoint center \
   -b 3000 -a 3000 \
   -R "$path_to_diffbind"diffenriched_h3k27ac_up.bed \
   -S "$path_to_bigwig"*.bw \
   --skipZeros -p 6 \
   -o "$path_to_project"matrix/matrix_diffenriched_h3k27ac_up_centered.gz \
   --outFileSortedRegions "$path_to_project"matrix/matrix_diffenriched_h3k27ac_centered_up.bed

   computeMatrix reference-point --referencePoint center \
   -b 3000 -a 3000 \
   -R "$path_to_diffbind"diffenriched_h3k27ac_down.bed \
   -S "$path_to_bigwig"*.bw \
   --skipZeros -p 6 \
   -o "$path_to_project"matrix/matrix_diffenriched_h3k27ac_down_centered.gz \
   --outFileSortedRegions "$path_to_project"matrix/matrix_diffenriched_h3k27ac_centered_down.bed

   computeMatrix reference-point --referencePoint center \
   -b 3000 -a 3000 \
   -R "$path_to_diffbind"diffenriched_h3k27ac_up_chr11.bed \
   -S "$path_to_bigwig"*.bw \
   --skipZeros -p 6 \
   -o "$path_to_project"matrix/matrix_diffenriched_h3k27ac_up_chr11_centered.gz \
   --outFileSortedRegions "$path_to_project"matrix/matrix_diffenriched_h3k27ac_centered_up_chr11.bed

   computeMatrix reference-point --referencePoint center \
   -b 3000 -a 3000 \
   -R "$path_to_diffbind"diffenriched_h3k27ac_down_chr11.bed \
   -S "$path_to_bigwig"*.bw \
   --skipZeros -p 6 \
   -o "$path_to_project"matrix/matrix_diffenriched_h3k27ac_down_chr11_centered.gz \
   --outFileSortedRegions "$path_to_project"matrix/matrix_diffenriched_h3k27ac_centered_down_chr11.bed

### Assigning H3K27Ac peaks to genes using Binding and Expression Target Analysis (BETA)
BETA basic -p "$path_to_chipr"GRANTA_H3K27Ac_chipr_all.bed --diff_expr "$path_to_deseq"genexpr_granta_names.txt --kind BSF -g hg38 --da 1000 -n basic \
            --gname2 -o "$path_to_project""beta/enhancers" --df 0.05 -c 0.5
                                          
### Assigning SEs to genes using Binding and Expression Target Analysis (BETA)
BETA basic -p "$path_to_rose_files"/rose_output_G2b/GRANTA_H3K27Ac_chipr_all_Gateway_SuperEnhancers.bed --diff_expr "$path_to_deseq"genexpr_granta_names.txt  \
            --kind BSF -g hg38 --da 1000 -n basic --gname2 -o "$path_to_project""beta/SEs" --df 0.05 -c 0.5
                                                           

