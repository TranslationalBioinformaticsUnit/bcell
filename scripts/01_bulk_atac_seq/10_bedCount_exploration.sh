#*************************************#
# ATAC-Seq pipeline - PAIR-END        #
# script: 10_bedCount_exploration.sh  #
# script: 10_bedCount_exploration.R   #
#*************************************#

# 10. Count the number of reads mapping to each feature - QUALITY CONTROL


#Input directory
In_Dir=$1
#Output directory
Out_Dir=$2
#Code directory
Code_Dir=$3
#Metadata file
meta_file=$4


# Exploration, normalization and filtering of count matrix.
/R/3.5.2/bin/R CMD BATCH --no-save --no-restore '--args filename="'$Out_Dir'/count_table.txt" filename2="'$meta_file'" workdir_bams="'$In_Dir'" workdir="'$Out_Dir'"' $Code_Dir/10_bedCount.R $Code_Dir/jobs/bedCount_output.out
