#*******************************#
# RNA-Seq pipeline - SINGLE-END #
#*******************************#


######################
##   CHECK INPUT    ##
######################
SamplesNames_Dir=$1 #File with names of samples to be analyzed
Fastq_Dir=$2 #Directory path where are the fastq files 
Work_Dir=$3 #Directory path where results are stored
Code_Dir=$4 #Directory where scripts are located
reference_genome=$5 #Reference genome 
Meta_File=$6 #Metadata file
Nodo=$7 #Define the node to run



######################
##     PIPELINE     ##
######################

# Run the pipeline for each fastq file
aux="hi"

for FastqID in $(cat $SamplesNames_Dir);
do
    echo "$FastqID"
 
    
    #########################################
    # Quality control of the original fastq #
    #########################################
    if [ ! -d $Work_Dir/01_fastQC ]; then
      mkdir $Work_Dir/01_fastQC
      chmod 777 $Work_Dir/01_fastQC
    fi
    #Parameters(4): FastqID Option(1) InputDir OutputDir
    ###Option=1 -> fastqc of original fastq
    ###Option=2 -> fastqc of trimmed fastq
    qsub -l h=$Nodo -hold_jid $aux -N QC_${FastqID} -e $Code_Dir/jobs/fastqc.stderr  -o $Code_Dir/jobs/fastqc.stdout $Code_Dir/01_fastqQC.sh $FastqID 1 $Fastq_Dir $Work_Dir/01_fastQC

	# MultiQC
	# when all the Fastqc were done, run the multiQC
	# multiqc $Work_Dir/01_fastQC/FastQC -o $Work_Dir/01_fastQC/MultiQC
 
 
    #################
    # Trimmed Reads #
    #################
    if [ ! -d $Work_Dir/02_trimmed_reads ]; then
      mkdir $Work_Dir/02_trimmed_reads
      chmod 777 $Work_Dir/02_trimmed_reads
    fi
    #Parameters(3): FastqID InputDir OutputDir
    qsub -l h=$Nodo -hold_jid QC_${FastqID} -N trim_${FastqID} -e $Code_Dir/jobs/trimmomatic.stderr  -o $Code_Dir/jobs/trimmomatic.stdout $Code_Dir/02_trimmedReads.sh $FastqID $Fastq_Dir $Work_Dir/02_trimmed_reads 
	
	
    ########################################
    # Quality control of the TRIMMED fastq #
    ########################################
    if [ ! -d $Work_Dir/02_trimmed_reads/trimmed_fastqQC ]; then
      mkdir $Work_Dir/02_trimmed_reads/trimmed_fastqQC
      chmod 777 $Work_Dir/02_trimmed_reads/trimmed_fastqQC
    fi
    #Parameters(4): FastqID Option(2) InputDir OutputDir
    ###Option=1 -> fastqc of original fastq
    ###Option=2 -> fastqc of trimmed fastq
    qsub -l h=$Nodo -hold_jid trim_${FastqID} -N QC_trim_${FastqID} -e $Code_Dir/jobs/fastqc.stderr  -o $Code_Dir/jobs/fastqc.stdout $Code_Dir/01_fastqQC.sh $FastqID 2 $Work_Dir/02_trimmed_reads $Work_Dir/02_trimmed_reads/trimmed_fastqQC

	# MultiQC
	# when all the trimmed fastqc were done, run the multiQC
	# multiqc $Work_Dir/02_trimmed_reads/trimmed_fastqQC -o $Work_Dir/02_trimmed_reads/trimmed_multiQC


    ################
    # Mapping STAR #
    ################
    if [ ! -d $Work_Dir/03_mapping_bam ]; then
      mkdir $Work_Dir/03_mapping_bam
      chmod 777 $Work_Dir/03_mapping_bam
    fi
    #Parameters(3): FastqID InputDir OutputDir
    qsub -l h=$Nodo -hold_jid QC_trim_${FastqID} -N mapping_${FastqID} -e $Code_Dir/jobs/bwamapping.stderr  -o $Code_Dir/jobs/bwamapping.stdout $Code_Dir/03.1_MappingSTAR.sh $FastqID $Work_Dir/02_trimmed_reads $Work_Dir/03_mapping_bam


    ##############################
    # Quality control of mapping #
    ##############################
    if [ ! -d $Work_Dir/03_mapping_bam/mappingQC ]; then
      mkdir $Work_Dir/03_mapping_bam/mappingQC
      chmod 777 $Work_Dir/03_mapping_bam/mappingQC
    fi
    #Parameters(4): FastqID OutputDir InputDir ReferenceGenome
    qsub -l h=$Nodo -hold_jid $aux -N mapQC_${FastqID} -e $Code_Dir/jobs/mappingQC.stderr  -o $Code_Dir/jobs/mappingQC.stdout $Code_Dir/03.2_MappingQC.sh $FastqID $Work_Dir/03_mapping_bam $Work_Dir/03_mapping_bam/mappingQC $reference_genome

	# MultiQC
	# when all the QC of mapping were done, run the multiQC
	# multiqc $Work_Dir/03_mapping_bam/mappingQC -o $Work_Dir/03_mapping_bam/MultiQC_bam


    ################
    # Count Matrix #
    ################
   if [ ! -d $Work_Dir/04_count_matrix ]; then
      mkdir $Work_Dir/04_count_matrix
      chmod 777 $Work_Dir/04_count_matrix
    fi   
    #Parameters(5): FastqID Mode_option InputDir OutputDir ReferenceDir
    ###Mode_option=1 -> intersection-nonempty(by default)
    ###Mode_option=2 -> union
    ###Mode_option=3 -> intersection_strict
    qsub -l h=$Nodo -N count_matrix_${FastqID} -e $Code_Dir/jobs/counts.stderr  -o $Code_Dir/jobs/counts.stdout $Code_Dir/04.1_count_table.sh $FastqID 1 $Work_Dir/03_mapping_bam $Work_Dir/04.1_count_matrix $reference_genome

aux='mapQC_'${FastqID}

done  


####################
# Merge count data #
####################
if [ ! -d $Work_Dir/04_count_matrix/merged_Count_tables]; then
	mkdir $Work_Dir/04_count_matrix/merged_Count_tables
    chmod 777 $Work_Dir/04_count_matrix/merged_Count_tables
fi   

#Parameters(4): In_Dir Out_Dir filename Code_Dir
qsub -l h=$Nodo -N merge_counts -e $Code_Dir/jobs/merge_counts.stderr  -o $Code_Dir/jobs/merge_counts.stdout $Code_Dir/04.2_merge_count_table.sh $Work_Dir/04_count_matrix $Work_Dir/04_count_matrix/merged_Count_tables Nonempty_merged_count_table $Code_Dir


############################
# Count Matrix Exploration #
############################
if [ ! -d $Work_Dir/04_count_matrix/exploration ]; then
	mkdir $Work_Dir/04_count_matrix/exploration
	chmod 777 $Work_Dir/04_count_matrix/exploration
fi   
#Exploration, normalization and filtering of count matrix.
#Parameters(4): In_Dir Out_Dir Code_Dir meta_file
qsub -l h=$Nodo -N counst_QC -e $Code_Dir/jobs/countsQC.stderr  -o $Code_Dir/jobs/countsQC.stdout $Code_Dir/05_count_table_exploration.sh $Work_Dir/04_count_matrix/exploration $Work_Dir/04_count_matrix/merged_Count_tables/Nonempty_merged_count_table.rda $Meta_File $Code_Dir 


###########################
# Differential Expression #
###########################
if [ ! -d $Work_Dir/05_differential_expression_analysis ]; then
	mkdir $Work_Dir/05_differential_expression_analysis
	chmod 777 $Work_Dir/05_differential_expression_analysis
fi   
#Differential analysis from normalized and filtered count (TMM) data.
#Parameters(4): Work_Dir rna_data meta_file Code_Dir
qsub -l h=$Nodo -N diff_expression -e $Code_Dir/jobs/diff_expression.stderr  -o $Code_Dir/jobs/diff_expression.stdout $Code_Dir/06_differential_expression_analysis.sh $Work_Dir/05_differential_expression_analysis $Work_Dir/04_count_matrix/exploration/voom_to_dea.RData $Meta_File $Code_Dir


#Differential expression of TFs
if [ ! -d $Work_Dir/05_differential_expression_analysis/tfs ]; then
	mkdir $Work_Dir/05_differential_expression_analysis/tfs
	chmod 777 $Work_Dir/05_differential_expression_analysis/tfs
fi
#Parameters(4): Work_Dir counts_file Code_Dir meta_file
qsub -l h=$Nodo -N diff_expression_tfs -e $Code_Dir/jobs/diff_expression_tfs.stderr  -o $Code_Dir/jobs/diff_expression_tfs.stdout $Code_Dir/07_differential_expression_tf.sh $Work_Dir/05_differential_expression_analysis/tfs $Work_Dir/04_count_matrix/exploration/voom_to_dea.RData $Meta_File $Code_Dir


##########################################################
# Transitional differencial analysis + Binary Clustering #
##########################################################
if [ ! -d $Work_Dir/06_transitional_differential_analysis ]; then
	mkdir $Work_Dir/06_transitional_differential_analysis
	chmod 777 $Work_Dir/06_transitional_differential_analysis
fi   
#Differential expression analysis through cell-stages transitions
#Parameters(4): Work_Dir rna_data meta_file Code_Dir
qsub -l h=$Nodo -N trans_diff_expression -e $Code_Dir/jobs/transitional_diff_expression.stderr  -o $Code_Dir/jobs/transitional_diff_expression.stdout $Code_Dir/8_transitional_differential_analysis.sh $Work_Dir/06_transitional_differential_analysis $Work_Dir/04_count_matrix/exploration/voom_to_dea.RData $Meta_File $Code_Dir


# Clustering of significant genes (binary data)
if [ ! -d $Work_Dir/06_transitional_differential_analysis/clustering ]; then
	mkdir $Work_Dir/06_transitional_differential_analysis/clustering
	chmod 777 $Work_Dir/06_transitional_differential_analysis/clustering
fi   
#Binary clustering from dea results.
#Parameters(3): sig_genes Work_Dir Code_Dir
qsub -l h=$Nodo -N trans_dea_clustering -e $Code_Dir/jobs/trans_dea_clustering.stderr  -o $Code_Dir/jobs/trans_dea_clustering.stdout $Code_Dir/09_transitional_deg_clustering.sh $Work_Dir/06_transitional_differential_analysis/dea_significant_binary.txt $Work_Dir/06_transitional_differential_analysis/clustering $Code_Dir


# Differential expression of TFs through cell-stages transitions
if [ ! -d $Work_Dir/06_transitional_differential_analysis/tfs ]; then
	mkdir $Work_Dir/06_transitional_differential_analysis/tfs
	chmod 777 $Work_Dir/06_transitional_differential_analysis/tfs
fi
#Parameters(5): Work_Dir rna_data meta_file tf_data Code_Dir
qsub -l h=$Nodo -N diff_expression_tfs -e $Code_Dir/jobs/diff_expression_tfs.stderr  -o $Code_Dir/jobs/diff_expression_tfs.stdout $Code_Dir/10_transitional_differential_analysis_tf.sh $Work_Dir/06_transitional_differential_analysis/tfs $Work_Dir/04_count_matrix/exploration/voom_to_dea.RData $Meta_File $Work_Dir/05_differential_expression_analysis/tfs/norm_filtered_tfs_data.txt $Code_Dir


# Clustering of significant TFs
if [ ! -d $Work_Dir/06_transitional_differential_analysis/tfs/clustering ]; then
	mkdir $Work_Dir/06_transitional_differential_analysis/tfs/clustering
	chmod 777 $Work_Dir/06_transitional_differential_analysis/tfs/clustering
fi   
#Binary clustering from dea results.
#Parameters(3): sig_genes Work_Dir Code_Dir
qsub -l h=$Nodo -N tf_clustering -e $Code_Dir/jobs/deg_clustering_tf.stderr  -o $Code_Dir/jobs/deg_clustering_tf.stdout $Code_Dir/11_transitional_deg_clustering_tf.sh $Work_Dir/06_transitional_differential_analysis/tfs/dea_significant_binary_tf.txt $Work_Dir/06_transitional_differential_analysis/tfs/clustering $Code_Dir    
