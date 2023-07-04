#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 06_FRiP_plot.sh           #
# script: 06_FRiP_plot.R            #
#***********************************#

# 6. Quality Control of Peak Calling - Fraction of Reads in Peaks (FRiP) PLOT


#Directory where scripts are located
Code_Dir=$1
#Input file with FRiP
Input_File=$2
#Work Directory
Work_Dir=$3


/R/3.5.2/bin/R CMD BATCH --no-save --no-restore '--args file="'${Input_File}'" workdir="'${Work_Dir}'"' $Code_Dir/06_FRiP_plot.R $Code_Dir/jobs/FRiP_output.out