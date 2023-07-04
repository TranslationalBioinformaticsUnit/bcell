#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 09_consensus_peaksQC.sh   #
# script: 09_jaccard_heatmap.R      #
#***********************************#

# 9. Quality Control of Consensus Peaks - Measuring dataset similarity


#Input directory BAM files
In_Dir=$1
#Output directory
Out_Dir=$2
#Code directory
Code_Dir=$3

# Measuring dataset similarity (Jaccard)
# https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md#bp6--measuring-dataset-similarity


file_labels=`ls $In_Dir | grep .bed`

file_labels=$file_labels!=consensus_cell_types.bed

    echo name" "$file_labels > $Out_Dir/pairwise_jaccard.txt
    for file1 in `ls $In_Dir | grep .bed`
    do
        echo -n $file1 >> $Out_Dir/pairwise_jaccard.txt

        for file2 in `ls $In_Dir | grep .bed`;
        do
            echo $file2
            # compute the jaccard stat for these two files.
            jaccard=`/bedtools2/bin/bedtools jaccard -a $In_Dir/$file1 -b $In_Dir/$file2`
            echo -n $jaccard
            # report the jaccard stat for these two files
            value_jaccard=$(echo $jaccard | cut -d " " -f 7)
            echo -n " "$value_jaccard >> $Out_Dir/pairwise_jaccard.txt
        done
        echo "\n" >> $Out_Dir/pairwise_jaccard.txt
    done

    # Heatmap of Jaccard
    /R/3.5.2/bin/R CMD BATCH --no-save --no-restore '--args filename="'$Out_Dir'/pairwise_jaccard.txt" workdir="'$Out_Dir'"' $Code_Dir/09_jaccard_heatmap.R $Code_Dir/jobs/consensus_output.out
