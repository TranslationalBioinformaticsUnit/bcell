#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 03_MappingQC.sh           #
#***********************************#

# 3. Mapping Quality Control


# Get sample name
FastqID=$1
#Output directory
Out_Dir=$2
#Input directory
In_Dir=$3


# Mark duplicates
java -jar -Xmx2g -Djava.io.tmpdir=`pwd`/tmp -jar /picard/build/libs/picard.jar \
		MarkDuplicates \
		INPUT=$In_Dir/${FastqID}.bam \
		OUTPUT=$Out_Dir/${FastqID}.marked.bam \
		METRICS_FILE=$Out_Dir/${FastqID}.marked.metrics.txt \
		REMOVE_DUPLICATES=false \
		ASSUME_SORTED=true \
		CREATE_INDEX=true

# Run samtools flagstats
/samtools-1.6/bin/samtools flagstat $Out_Dir/${FastqID}.marked.bam > $Out_Dir/${FastqID}.marked.flagstat.txt
