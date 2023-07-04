#***************************************#
# RNA-Seq pipeline - SINGLE-END         #
# script: 06_differential_expression.sh #
# script: 06_differential_expression.R  #
#***************************************#

# 6. Differential Expression Analysis


#Work directory
Work_Dir=$1
#Count matrix
rna_data=$2
#metadata
metadata=$3
#Code directory
Code_Dir=$4


#Run differential expression analysis
/R/3.5.2/bin/R CMD BATCH --no-save --no-restore '--args workdir="'$Work_Dir'" rna_data="'$rna_data'" metadata="'$metadata'"' $Code_Dir/06_differential_expression.R $Code_Dir/jobs/differential_expression.out
