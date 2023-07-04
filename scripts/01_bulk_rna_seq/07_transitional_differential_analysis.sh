#**************************************************#
# RNA-Seq pipeline - SINGLE-END                    #
# script: 07_transitional_differential_analysis.sh #
# script: 07_transitional_differential_analysis.R  #
#**************************************************#

# 8. Differential Expression Analysis for cell transitions


#Work directory
Work_Dir=$1
#Count matrix
rna_data=$2
#metadata
metadata=$3
#Code directory
Code_Dir=$4


# Transitional Differential Analysis
/R/3.5.2/bin/R CMD BATCH --no-save --no-restore '--args workdir="'$Work_Dir'" rna_data="'$rna_data'" metadata="'$metadata'"' $Code_Dir/07_transitional_differential_analysis.R $Code_Dir/jobs/transitional_differential_analysis.out
