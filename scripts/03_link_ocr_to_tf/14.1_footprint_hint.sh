#************************************#
# FOOTPRINT IDENTIFICATION           #
# script: 14.1_footprint_hint.sh     #
# script: 14.1_histogram_footprint.R #
#************************************#


# 14.1 Footprint - HINT (http://www.regulatory-genomics.org/hint/introduction/)


#File_Input (BAM)
Input_BAM=$1
#File_Input (BAM)
Input_BED=$2
#Ref_genom
genome=$3
#Output directory
Out_Dir=$4
#Output name
Sample_Name=$5
#Meta_File
Meta_file=$6



#Run RGT-HINT 
/python/bin/rgt-hint footprinting --atac-seq --paired-end --organism=$genome --output-prefix=${Sample_Name}_footprint --output-location=$Out_Dir $Input_BAM $Input_BED

#Filter noisy footprints with low counts. Threshold here is 20 counts. You can set to your own threshold. Maybe make a histogram of the footprint distributions.
#In our sample test, 1Q is arround 50 in median. So let the threshold at 20 counts.

# Inspect the 1Q of all samples with R 
/R/3.5.2/bin/R CMD BATCH --no-save --no-restore '--args filename="'$Meta_file'" workdir="'$Out_Dir'"' $Code_Dir/14.1_histogram_footprint.R $Code_Dir/jobs/histogram_footprint_output.out

#Set a filter to 20 counts per footprint
awk -v OFS='\t' '$5>=20{print $1,$2,$3,$4,$5,$6}' $Out_Dir/${Sample_Name}_footprint.bed > $Out_Dir/${Sample_Name}_footprint_20.bed
