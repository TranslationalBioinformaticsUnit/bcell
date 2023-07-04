#****************************************************#
# ATAC-Seq pipeline - PAIR-END                       #
# script: 13_diff_binding_transitions_clustering.sh  #
# script: 13_diff_binding_transitions_clustering.R   #
#****************************************************#

# 13. Peak clustering


#Significant_peaks
sig_peaks=$1
#Peak annotation file
peak_anno=$2
#Output directory
Work_Dir=$3
#Code directory
Code_Dir=$4


# Clustering
/R/3.5.2/bin/R CMD BATCH --no-save --no-restore '--args sig_peaks="'$sig_peaks'" peaks_gene_file="'$peak_anno'"workdir="'$Work_Dir'"' $Code_Dir/13_diff_binding_transitions_clustering.R $Code_Dir/jobs/diff_binding_clustering.out
