#********************************#
# FOOTPRINT TO TF MOTIF          #
# script: 15_vierstra_motifs.sh  #
#********************************#


# 15. Get the candidate motifs on our OCRs based on Vierstra database (https://www.vierstra.org/)

# We first need to format our consensus OCRs BED file to give as input for tabix (ex. chr4:144021505-144021837)


#File_Input (BED)
Peaks_file=$1
#Work_directory
WorkDir=$2


cd $WorkDir


##Get the candidate motifs on peak regions
for peak in $(cat $Peaks_file);
do
    echo "$peak"
    chr=$(echo $peak | cut -d ";" -f1)
    echo "$chr"

    tabix https://resources.altius.org/~jvierstra/projects/motif-clustering/releases/v1.0/hg38.archetype_motifs.v1.0.bed.gz ${chr} >> ${WorkDir}/reference_vierstra_motifs.txt 

done  
