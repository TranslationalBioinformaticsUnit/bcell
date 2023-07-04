#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 03_MappingBWA.sh          #
#***********************************#

# 3. Map reads to reference genome using BWA


#Get sample name
FastqID=$1
#Output directory
Out_Dir=$2
#Input directory
In_Dir=$3
#Get the option
Option=$4
#Reference genome
Ref_genome=$5 


if (( $Option == 1 )) 
then
# Map paired-end reads and add read group info.
/bwa/bwa mem -t 20 $Ref_genome  \
                   $In_Dir/${FastqID}_1_trimmed.fastq.gz \
                   $In_Dir/${FastqID}_2_trimmed.fastq.gz \
                   -R '@RG\tID:'${FastqID}'\tSM:'${FastqID}'\tPL:Illumina\tCN:CBRA\tLB:Fragment' \
                   | /samtools sort -@20 -o $Out_Dir/${FastqID}.bam - 
 
else
# Map single-end reads and add read group info
/bwa/bwa mem -t 10 $Ref_genome \
                   $In_Dir/${FastqID}_1_trimmed.fastq.gz \
                   -R '@RG\tID:'${FastqID}'\tSM:'${FastqID}'\tPL:Illumina\tCN:CBRA\tLB:Fragment\tPU:'${FastqID} \
                   | /opt_repository/samtools sort -@10 -o $Out_Dir/${FastqID}.bam -                        
fi

# Index BAM file
/samtools-1.6/bin/samtools index $Out_Dir/${FastqID}.bam

# Run samtools flagstats
/samtools-1.6/bin/samtools flagstat $Out_Dir/${FastqID}.bam > $Out_Dir/${FastqID}.flagstat.txt
