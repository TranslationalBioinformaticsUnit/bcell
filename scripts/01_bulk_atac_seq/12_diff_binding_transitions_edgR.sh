#**********************************************#
# ATAC-Seq pipeline - PAIR-END                 #
# script: 12_diff_binding_transitions_edgR.sh  #
# script: 12_diff_binding_transitions_edgR.R   #
#**********************************************#

# 12. Differential Binding Analysis


#Work directory
Work_Dir=$1
#Counts file
filename=$2
#Code directory
Code_Dir=$3
#Metadata file
meta_file=$4


# Differential Binding Analysis through transitions
/R/3.5.2/bin/R CMD BATCH --no-save --no-restore '--args filename="'$filename'" meta_file="'$meta_file'" workdir="'$Work_Dir'"' $Code_Dir/12_diff_binding_transitions_edgR.R $Code_Dir/jobs/diffBind_transitions_output.out
