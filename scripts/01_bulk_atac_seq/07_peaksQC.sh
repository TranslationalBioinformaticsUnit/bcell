#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 07_peaksQC.sh             #
# script: 07_jaccard_heatmap.R      #
#***********************************#

# 7. Quality Control of Peak Calling - Measuring dataset similarity


#File with samples (metadata)
In_File=$1
#Output directory
Out_Dir=$2
#Code directory
Code_Dir=$3
#Output file
Out_File=$4


declare arrays
declare -a file_names
declare -a sample_names
i=1

for BED_file in $(cat $In_File);
do

  file_names[$i]=$(echo $BED_file | cut -d ";" -f6)
  sample_names[$i]=$(echo $BED_file | cut -d ";" -f1)

  i=$((i+1))
done

# Measuring dataset similarity (Jaccard)
# https://github.com/arq5x/bedtools-protocols/blob/master/bedtools.md#bp6--measuring-dataset-similarity

   echo ${sample_names[*]} > $Out_Dir/$Out_File

#Remove header names
unset file_names[1]
unset sample_names[1]
    
    i=2
    for file1 in "${file_names[@]}";
    do
    echo -n ${sample_names[$i]} >> $Out_Dir/$Out_File
    
        for file2 in "${file_names[@]}";
        do
        echo $file2
            # compute the jaccard stat for these two files.
            jaccard=`/bedtools2/bin/bedtools jaccard -a $file1 -b $file2`
       
            # report the jaccard stat for these two files
            value_jaccard=$(echo $jaccard | cut -d " " -f 7)
            echo -n " "$value_jaccard >> $Out_Dir/$Out_File

        done
            
    i=$((i+1))
    echo "" >> $Out_Dir/$Out_File
    done 

    
    # Heatmap of Jaccard
    /R/3.5.2/bin/R CMD BATCH --no-save --no-restore '--args filename="'${Out_Dir}'/'${Out_File}'" filename2="'${In_File}'" workdir="'${Out_Dir}'"' $Code_Dir/07_jaccard_heatmap.R $Code_Dir/jobs/jaccard_output.out
