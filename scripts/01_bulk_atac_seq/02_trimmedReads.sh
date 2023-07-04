#*******************************#
# ATAC-Seq pipeline - PAIR-END  #
# script: 02_trimmedReads.sh    #
#*******************************#

# 2. Trim reads and remove adapters with Trimmomatic


# Get sample name
FastqID=$1
#Get output directory
Out_Dir=$2
#Get input directory
In_Dir=$3
# Get the option
Option=$4


if (( $Option == 1 )) 
then
# Run Trimmomatic for paired-end reads (PE)
java -jar trimmomatic/classes/trimmomatic.jar PE -threads 20 -phred33 $In_Dir/${FastqID}_R1_001.fastq.gz $In_Dir/${FastqID}_R2_001.fastq.gz $Out_Dir/${FastqID}_1_trimmed.fastq.gz $Out_Dir/${FastqID}_U1_trimmed.fastq.gz $Out_Dir/${FastqID}_2_trimmed.fastq.gz $Out_Dir/${FastqID}_U2_trimmed.fastq.gz ILLUMINACLIP:/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15

else
# Run Trimmomatic for single-end reads (SE)
java -jar /trimmomatic/classes/trimmomatic.jar SE -threads 20 -phred33 $In_Dir/${FastqID}_R1_001.fastq.gz $Out_Dir/${FastqID}_1_trimmed.fastq.gz ILLUMINACLIP:/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:15
fi
