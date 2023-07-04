#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 11_peak_annotation.sh     #
# script: 11_peak_annotation.R      #
#***********************************#

# 11. Genomic mapping/annotation of OCRs


#Input BED file
Input_BED=$1
#Work directory
Work_Dir=$2
#Count matrix
Count_table=$3
#Code directory
Code_Dir=$4


/R/3.5.2/bin/R CMD BATCH --no-save --no-restore '--args workdir="'$Work_Dir'" peaks_file="'$Input_BED'" count_matrix="'$Count_table'"' $Code_Dir/11_peak_annotation.R $Code_Dir/jobs/annotation_output.out
