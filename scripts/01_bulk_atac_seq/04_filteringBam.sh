#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 04_filteringBam.sh        #
#***********************************#

# 4. Filtering BAM files


# Get sample name
FastqID=$1
#Output directory
Out_Dir=$2
#Input directory
In_Dir=$3
#Blacklist
BLACKLIST=$4


# Run some BAM filtering (pcr duplicates, reads unmapped, mapq quality below 10 and blacklist)
/samtools-1.6/bin/samtools view -b -F 4 -F 1024 -q 10 $In_Dir/${FastqID}.marked.bam | /samtools-1.6/bin/samtools view -b -o $Out_Dir/${FastqID}.blacklist.bam -U $Out_Dir/${FastqID}.filtered.bam -L $BLACKLIST

# Index filtered bam file
/samtools-1.6/bin/samtools index $Out_Dir/${FastqID}.filtered.bam

# Keep just known chromosomes, remove random, unknown and chrM
/samtools-1.6/bin/samtools view -b $Out_Dir/${FastqID}.filtered.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY > $Out_Dir/${FastqID}.clean.bam

# Index filtered bam file
/samtools-1.6/bin/samtools index $Out_Dir/${FastqID}.clean.bam

# Picard statistics
java -jar -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar /picard/build/libs/picard.jar \
		MarkDuplicates \
		INPUT=$Out_Dir/${FastqID}.clean.bam \
		OUTPUT=$Out_Dir/${FastqID}.clean.marked.bam \
		METRICS_FILE=$Out_Dir/${FastqID}.clean.marked.metrics.txt

# Run samtools flagstats
/samtools-1.6/bin/samtools flagstat $Out_Dir/${FastqID}.clean.marked.bam > $Out_Dir/${FastqID}.clean.marked.flagstat.txt