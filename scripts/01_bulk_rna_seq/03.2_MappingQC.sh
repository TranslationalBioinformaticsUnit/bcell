#*******************************#
# RNA-Seq pipeline - SINGLE-END #
# script: 03.2_MappingQC.sh     #
#*******************************#

# 3.2. Quality Control of STAR mapping


# Get the fastq id file
FastqID=$1
#Get input directory
Input_Dir=$2
#Get output directory
Output_Dir=$3
#Get reference directory
Reference_Dir=$4

/java -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar /picard/build/libs/picard.jar CollectAlignmentSummaryMetrics \
		I= $Input_Dir/${FastqID}.dedup.bam \
		O= $Output_Dir/${FastqID}.metrics.txt \
		REFERENCE_SEQUENCE= /Reference_genome_Star/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
		VALIDATION_STRINGENCY=LENIENT
		