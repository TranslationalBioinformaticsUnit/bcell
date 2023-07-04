#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 06_FRiP_calculation.sh    #
#***********************************#

# 6. Quality Control of Peak Calling - Fraction of Reads in Peaks (FRiP)


# Get sample name
FastqID=$1
#Input directory BAM files
In_Dir=$2
#Input directory BED files
In_Dir2=$3
#Output directory
Out_Dir=$4


# FRiP calculation
total=$(/samtools-1.6/bin/samtools view -c $In_Dir/${FastqID}.clean.bam)

reads_in_peaks=$(/samtools-1.6/bin/samtools view -c $In_Dir/${FastqID}.clean.bam -L $In_Dir2/${FastqID}.peak_calling_peaks.narrowPeak)

reads_in_filtered_peaks=$(/samtools-1.6/bin/samtools view -c $In_Dir/${FastqID}.clean.bam -L $In_Dir2/${FastqID}.peak_calling.narrowPeak.filtered.bed)

peaks=$(wc -l $In_Dir2/${FastqID}.peak_calling.narrowPeak.filtered.bed)
peaks_filtered=$(wc -l $In_Dir2/${FastqID}.peak_calling_peaks.narrowPeak)

echo "$FastqID $total $reads_in_peaks $reads_in_filtered_peaks $peaks $peaks_filtered" >> $Out_Dir/frip_input.txt
