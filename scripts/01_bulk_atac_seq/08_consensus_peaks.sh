#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 08_consensus_peaks.sh     #
# script: 08_consensus_peaks.R      #
#***********************************#

# 8. Consensus peaks for each cell type using GenomicRanges and rtracklayer
#https://ro-che.info/articles/2018-07-11-chip-seq-consensus


#Code directory
Code_Dir=$1
#File input
In_File=$2
#Work directory
Work_Dir=$3


/R/3.5.2/bin/R CMD BATCH --no-save --no-restore '--args filename="'${In_File}'" filestats="'${Work_Dir}'/04_merged_filtered_bam/multiQC/multiqc_data/multiqc_general_stats.txt" workdir="'${Work_Dir}'/08_consensus_peaks"' $Code_Dir/08_consensus_peaks.R $Code_Dir/jobs/consensus_output.out
