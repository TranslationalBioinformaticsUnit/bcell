#*******************************#
# ATAC-Seq pipeline - PAIR-END  #
# script: 01_fastqQC.sh			#
#*******************************#

# 1. Quality control of a fastq file


# Get sample name
FastqID=$1
# Get the option
Option=$2
#Get output directory
Out_Dir=$3
#Get input directory
In_Dir=$4


if (( $Option == 1 )) 
then
  fastqc -o $Out_Dir/ -t 10 $In_Dir/${FastqID}_R1_001.fastq.gz $In_Dir/${FastqID}_R2_001.fastq.gz
fi

if (( $Option == 2 ))
then 
  fastqc -o $Out_Dir/ -t 10 $In_Dir/${FastqID}_1_trimmed.fastq.gz $In_Dir/${FastqID}_2_trimmed.fastq.gz
  fastqc -o $Out_Dir/ -t 10 $In_Dir/${FastqID}_U1_trimmed.fastq.gz $In_Dir/${FastqID}_U2_trimmed.fastq.gz
fi
