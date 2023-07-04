#********************************************#
# FOOTPRINT IDENTIFICATION                   #
# script: 14.3_unify_footprint.sh            #
# script: 14.3_unify_footprint.R             #
# script: 14.4_consensus_unify_footprint.R   #
#********************************************#


# 14.3 Merge footprints from HINT and Wellington outputs and get the consensus footprint for each cell type


#Work directory
Out_Dir=$1
#Code directory
Code_Dir=$2
#Meta_File
Meta_file=$3


#Get the unified footprint (hint + wellington approaches) for each sample
/R/3.5.2/bin/R CMD BATCH --no-save --no-restore '--args filename="'$Meta_file'" workdir="'$Out_Dir'"' $Code_Dir/14.3_unify_footprint.R $Code_Dir/jobs/unify_footprint_output.out

#Get the consensus footprint for each cell type
/R/3.5.2/bin/R CMD BATCH --no-save --no-restore '--args filename="'$Meta_file'" workdir="'$Out_Dir'"' $Code_Dir/14.4_consensus_unified_footprint.R $Code_Dir/jobs/consensus_footprint_output.out
