#*******************************************#
# RNA-Seq pipeline - SINGLE-END             #
# script: 08_transitional_deg_clustering.sh #
# script: 08_transitional_deg_clustering.R  #
#*******************************************#

# 9. Gene clustering


#transition_Significant_genes
sig_genes=$1
#Output directory
Work_Dir=$2
#Code directory
Code_Dir=$3


# Clustering
/R/3.5.2/bin/R CMD BATCH --no-save --no-restore '--args trans_dea_sig="'$sig_genes'" workdir="'$Work_Dir'"' $Code_Dir/08_transitional_deg_clustering.R $Code_Dir/jobs/transitional_deg_clustering.out
