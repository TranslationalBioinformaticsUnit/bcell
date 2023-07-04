#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 05_macs2.sh               #
#***********************************#

# 5. PeakCalling for each samples using MACS2


# Get sample name
FastqID=$1
#Output directory
Out_Dir=$2
#Input directory
In_Dir=$3
#Option
Option=$4


if (( $Option == 1 )) 
then

/python/bin/macs2 callpeak -t $In_Dir/${FastqID}.clean.bam --outdir $Out_Dir -g hs -f BAMPE -q 0.01 --nomodel --shift -100 --extsize 200 -B --SPMR -n ${FastqID}.peak_calling

else

/python/bin/macs2 callpeak -t $In_Dir/${FastqID}.clean.bam --outdir $Out_Dir -g hs -f BAM -q 0.01 --nomodel --shift -100 --extsize 200 -B --SPMR -n ${FastqID}.peak_calling

fi


#Filtering those peaks out of a logical range
awk '{if ($3-$2 > 20) print $0;}' $Out_Dir/${FastqID}.peak_calling_peaks.narrowPeak > $Out_Dir/${FastqID}.peak_calling.narrowPeak.filtered.bed
