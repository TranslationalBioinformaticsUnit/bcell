#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 04_mergeBam.sh            #
#***********************************#

# 4. Merge BAMs


SamplesNames_Dir=$1 #File with paired names of samples to be merged
Work_Dir=$2 #Directory to work
Code_Dir=$3 #Directory where scripts are located
Bam_Dir=$4 #Directory were BAM files are located
Nodo=$5 #Define the node to run


if [ ! -d $Work_Dir/04_merged_filtered_bam ]; then
  mkdir $Work_Dir/04_merged_filtered_bam
  chmod 777 $Work_Dir/04_merged_filtered_bam
fi  


aux="hi"

for FastqID in $(cat $SamplesNames_Dir);
do
    FastqID_1=$(echo $FastqID | cut -d ";" -f1)
    echo "$FastqID_1"
    FastqID_2=$(echo $FastqID | cut -d ";" -f2)
    echo "$FastqID_2"

    #Merge Bam
    #Parameters(4): FastqID_1 FastqID_2
    qsub -l h=$Nodo -hold_jid $aux -N mergeBam_${FastqID_1} -e $Code_Dir/jobs/mergeBam.stderr  -o $Code_Dir/jobs/mergeBam.stdout $Code_Dir/04_mergeBam_internal.sh ${FastqID_1} ${FastqID_2} $Bam_Dir $Work_Dir/04_merged_filtered_bam
  
    aux='mergeBam_'${FastqID_1}   
done
