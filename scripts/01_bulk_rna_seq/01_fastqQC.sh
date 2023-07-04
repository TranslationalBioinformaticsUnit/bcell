#*******************************#
# RNA-Seq pipeline - SINGLE-END #
# script: 01_fastQC.sh			#
#*******************************#

# 1. Quality control of a fastq file


# Get the fastq id file
FastqID=$1
# Get the option
Option=$2
#Get input direction
Input_Dir=$3
#Get output direction
Output_Dir=$4

if (( $Option == 1 )) 
then
  # Run quality control of the original fastq using FastQC
  fastqc -o $Output_Dir/ -t 10 $Input_Dir/${FastqID}.fastq.gz
else
  #Quality Control del trimmed fastq
  fastqc -o $Output_Dir/ -t 10 $Input_Dir/${FastqID}_trimmed.fastq.gz
fi