#*******************************#
# RNA-Seq pipeline - SINGLE-END #
# script: 04.1_count_table.sh   #
#*******************************#

# 4.1. Count the number of reads mapping to each feature
# https://htseq.readthedocs.io/en/release_0.11.1/count.html


# Get the fastq id file
FastqID=$1
#Option
Option=$2
#Get input directory
Input_Dir=$3
#Get output directory
Output_Dir=$4
#Get reference directory
Reference_genome=$5


#htseq-count
if (( $Option == 1 )) 
then
  /python/bin/htseq-count -i gene_id --mode=intersection-nonempty --nonunique=none --format=bam $Input_Dir/${FastqID}.dedup.bam $Reference_genome> $Output_Dir/${FastqID}_nonempty_count_table.txt
elif (( $Option == 2 )) 
then
  /python/bin/htseq-count -i gene_id --mode=union --nonunique=none --format=bam $Input_Dir/${FastqID}.dedup.bam $Reference_genome > $Output_Dir/${FastqID}_union_count_table.txt
elif (( $Option == 3 )) 
then
  /python/bin/htseq-count -i gene_id --mode=intersection-strict --nonunique=none --format=bam $Input_Dir/${FastqID}.dedup.bam $Reference_genome > $Output_Dir/${FastqID}_intersection-strict_count_table.txt
fi
