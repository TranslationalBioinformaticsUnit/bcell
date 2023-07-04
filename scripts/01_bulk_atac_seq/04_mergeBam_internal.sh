#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 04_mergeBam_internal.sh   #
#***********************************#

# 4. Merge BAMs files 


# Get sample name 1
FastqID_1=$1
# Get sample name 1
FastqID_2=$2
# Output directory
In_Dir=$3
# Output directory
Out_Dir=$4


# Merge samples and store in 04_merged_filtered_bam 
/samtools-1.6/bin/samtools merge $Out_Dir/${FastqID_1}.clean.bam $In_Dir/${FastqID_1}.clean.bam $In_Dir/${FastqID_2}.clean.bam

# Index filtered bam file
/samtools-1.6/bin/samtools index $Out_Dir/${FastqID_1}.clean.bam

# Run samtools flagstats
/samtools-1.6/bin/samtools flagstat $Out_Dir/${FastqID_1}.clean.bam > $Out_Dir/${FastqID_1}.merge.flagstat.txt
