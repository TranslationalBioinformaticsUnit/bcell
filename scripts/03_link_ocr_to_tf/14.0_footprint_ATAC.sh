#*******************************#
# ATAC-Seq - FOOTPRINT          #
#*******************************#


######################
##   CHECK INPUT    ##
######################
SamplesNames_Dir=$1 #File with names of samples to be analyzed
Fastq_Dir=$2 #Directory path where are the fastq files 
Work_Dir=$3 #Directory path where results are stored
Code_Dir=$4 #Directory where scripts are located
Ref_genome=$5 #Reference genome file
GTF_genome=$6 #GTF genome for annotation
BlackList=$7 #Blacklist BED file 
Meta_File=$8 #Metadata file
Nodo=$9 #Define the node to run



################################
##  FOOTPRINT IDENTIFICATION  ##
################################


###############################
# Footprinting HINT algorithm #
###############################
if [ ! -d $Work_Dir/12_footprinting ]; then
	mkdir $Work_Dir/12_footprinting
    chmod 777 $Work_Dir/12_footprinting
fi

if [ ! -d $Work_Dir/12_footprinting/hint ]; then
	mkdir $Work_Dir/12_footprinting/hint
    chmod 777 $Work_Dir/12_footprinting/hint
fi

aux="hi"

for FastqID in $(cat $Meta_File);
do
    bam_file=$(echo $FastqID | cut -d ";" -f5)
    echo "$bam_file"
    sample_name=$(echo $FastqID | cut -d ";" -f1)
    echo "$sample_name"
    bed_file=$(echo $FastqID | cut -d ";" -f6)
    echo "$bed_file"

    qsub -l h=$Nodo -hold_jid $aux -N footprinting_${sample_name} -e $Code_Dir/jobs/footprint.stderr  -o $Code_Dir/jobs/footprint.stdout $Code_Dir/14.1_footprint_hint.sh $bam_file $bed_file hg38 $Work_Dir/12_footprinting/hint $sample_name $Meta_File

    aux='footprinting_'${sample_name}
done



#####################################
# Footprinting Wellington algorithm #
#####################################
if [ ! -d $Work_Dir/12_footprinting/wellington ]; then
	mkdir $Work_Dir/12_footprinting/wellington
	chmod 777 $Work_Dir/12_footprinting/wellington
fi

aux="hi"

for FastqID in $(cat $Meta_File);
do
    bam_file=$(echo $FastqID | cut -d ";" -f5)
    echo "$bam_file"
    sample_name=$(echo $FastqID | cut -d ";" -f1)
    echo "$sample_name"
    bed_file=$(echo $FastqID | cut -d ";" -f8)
    echo "$bed_file"

    if [ ! -d $Work_Dir/13_footprint_wellington/$sample_name ]; then
      mkdir $Work_Dir/13_footprint_wellington/$sample_name
      chmod 777 $Work_Dir/13_footprint_wellington/$sample_name
    fi

    qsub -l h=$Nodo -hold_jid $aux -N footprinting2_${sample_name} -e $Code_Dir/jobs/footprint_wellington.stderr  -o $Code_Dir/jobs/footprint_wellington.stdout $Code_Dir/14.2_footprint_wellington.sh $bam_file $bed_file $Work_Dir/12_footprinting/wellington

    aux='footprinting2_'${sample_name}
done


#########################################
# Footprinting UNIFY: HINT + Wellington #
#########################################
if [ ! -d $Work_Dir/12_footprinting/unified ]; then
	mkdir $Work_Dir/12_footprinting/unified
	chmod 777 $Work_Dir/12_footprinting/unified
fi

qsub -l h=$Nodo -N unify_footprint -e $Code_Dir/jobs/footprint_unify.stderr  -o $Code_Dir/jobs/footprint_unify.stdout $Code_Dir/14.3_unify_footprint.sh $Work_Dir/12_footprinting/unified $Code_Dir $Meta_File 


