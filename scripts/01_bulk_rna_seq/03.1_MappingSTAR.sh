#*******************************#
# RNA-Seq pipeline - SINGLE-END #
# script: 03.1_MappingSTAR.sh   #
#*******************************#

# 3.1. Mapping with STAR program


#Fastq ID
FASTQ_ID=$1
#Directory with fastaq files
FASTQ_DIR=$2
#Outpur directory
BAM_DIR=$3


# Reference dir: STAR_INDEX
STAR_INDEX=/Reference_genome_Star/


/STAR-2.6.1c/source/STAR --runThreadN 8 --genomeDir $STAR_INDEX --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFileNamePrefix "${BAM_DIR}/${FASTQ_ID}." --readFilesIn "${FASTQ_DIR}/${FASTQ_ID}_trimmed.fastq.gz"

/samtools-1.6/bin/samtools index "${BAM_DIR}/${FASTQ_ID}.Aligned.sortedByCoord.out.bam"

/java -jar -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar /picard/build/libs/picard.jar MarkDuplicates I="${BAM_DIR}/${FASTQ_ID}.Aligned.sortedByCoord.out.bam" O="${BAM_DIR}/${FASTQ_ID}.dedup.bam" M="${BAM_DIR}/${FASTQ_ID}.marked_dup_metrics.txt";

/samtools-1.6/bin/samtools index "${BAM_DIR}/${FASTQ_ID}.dedup.bam";
