#***********************************#
# RNA-Seq pipeline - SINGLE-END     #
# script: 04.2_merge_count_table.sh #
# script: 04.2_merge_count_table.R  #
#***********************************#

# 4.2. Generate the count matrix with all the samples
# https://htseq.readthedocs.io/en/release_0.11.1/count.html


#Input directory
In_Dir=$1
#Output directory
Out_Dir=$2
#Output matrix File name
filename=$3
#Code directory
Code_Dir=$4

# Compute count matrix 
/R/3.5.2/bin/R CMD BATCH --no-save --no-restore '--args indir="'$In_Dir'" outdir="'$Out_Dir'" filename="'$filename'"' $Code_Dir/04.2_merge_count_table.R $Code_Dir/jobs/merge_count_table.out
