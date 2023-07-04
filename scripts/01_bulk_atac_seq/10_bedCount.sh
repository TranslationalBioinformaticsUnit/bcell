#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 10_bedCount.sh            #
#***********************************#

# 10. Count the number of reads mapping to each feature (OCR)


#Input directory
In_Dir=$1
#Input directory
bed_file=$2
#Output directory
Out_Dir=$3


filenames=$(ls -d $In_Dir/{*,.*}| grep .clean.bam$)

/bedtools2/bin/bedtools multicov -bams $filenames -bed $bed_file -D > $Out_Dir/count_table.txt
