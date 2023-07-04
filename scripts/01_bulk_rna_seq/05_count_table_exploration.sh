#***************************************#
# RNA-Seq pipeline - SINGLE-END         #
# script: 05_count_table_exploration.sh #
# script: 05_count_table_exploration.R  #
#***************************************#

# 5. Count matrix exploration


#Work directory
Work_Dir=$1
#Count matrix
count_matrix=$2
#metadata
metadata=$3
#Code directory
Code_Dir=$4


# Compute count matrix 
/R/3.5.2/bin/R CMD BATCH --no-save --no-restore '--args workdir="'$Work_Dir'" count_table="'$count_matrix'" metadata="'$metadata'"' $Code_Dir/05_count_table_exploration.R $Code_Dir/jobs/count_table_exploration.out
