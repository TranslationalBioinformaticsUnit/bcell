#*******************************#
# RNA-Seq pipeline - SINGLE-END #
# script: 02_trimmedRead.sh     #
#*******************************#

# 2. Trimming of a fastq file


# Get the fastq id file
FastqID=$1
#Get input direction
Input_Dir=$2
#Get output direction
Output_Dir=$3

/java -jar /trimmomatic/classes/trimmomatic.jar SE -threads 20 -phred33 $Input_Dir/${FastqID}.fastq.gz $Output_Dir/${FastqID}_trimmed.fastq.gz ILLUMINACLIP:/trimmomatic/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:30
