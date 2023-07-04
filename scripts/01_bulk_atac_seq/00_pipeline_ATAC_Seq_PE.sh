#*******************************#
# ATAC-Seq pipeline - PAIR-END  #
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



######################
##     PIPELINE     ##
######################

## Run the pipeline for each fastq file
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
    #Parameters(4): FastqID Option(1) OutputDir InputDir
    ###Option=1 -> fastqc of original fastq
    ###Option=2 -> fastqc of trimmed fastq
    qsub -l h=$Nodo -hold_jid $aux -N QC_${FastqID} -e $Code_Dir/jobs/fastqc.stderr  -o $Code_Dir/jobs/fastqc.stdout $Code_Dir/01_fastqQC.sh $FastqID 1 $Work_Dir/01_fastQC $Fastq_Dir

 
    #################
    # Trimmed Reads #
    #################
    if [ ! -d $Work_Dir/02_trimmed_reads ]; then
      mkdir $Work_Dir/02_trimmed_reads
      chmod 777 $Work_Dir/02_trimmed_reads
#   fi
    #Parameters(4): FastqID InputDir OutputDir Option(1)
    ###Option=1 -> fastqc of PE
    ###Option=2 -> fastqc of SE
     qsub -l h=$Nodo -hold_jid QC_${FastqID} -N trim_${FastqID} -e $Code_Dir/jobs/trimmomatic.stderr  -o $Code_Dir/jobs/trimmomatic.stdout $Code_Dir/02_trimmedReads.sh $FastqID $Work_Dir/02_trimmed_reads $Fastq_Dir 1

        ##############
        # ALTERNATIVE: Trimmed Reads with cut reads (CROP defined at 38) [for second batch reads]
        ##############
#       if [ ! -d $Work_Dir/02_trimmed_reads ]; then
#        mkdir $Work_Dir/02_trimmed_reads
#        chmod 777 $Work_Dir/02_trimmed_reads
#      	fi
		#Parameters(4): FastqID InputDir OutputDir Option(1)
		###Option=1 -> fastqc of PE
		###Option=2 -> fastqc of SE
#       qsub -l h=$Nodo -hold_jid QC_${FastqID} -N trim_${FastqID} -e $Code_Dir/jobs/trimmomatic.stderr  -o $Code_Dir/jobs/trimmomatic.stdout $Code_Dir/02_trimmedReads_cut.sh $FastqID $Work_Dir/02_trimmed_reads $Fastq_Dir 1

        
    ########################################
    # Quality control of the TRIMMED fastq #
    ########################################
    if [ ! -d $Work_Dir/02_trimmed_reads/trimQC ]; then
      mkdir $Work_Dir/02_trimmed_reads/trimQC
      chmod 777 $Work_Dir/02_trimmed_reads/trimQC
    fi
    #Parameters(4): FastqID Option(2) OutputDir InputDir
    ###Option=1 -> fastqc of original fastq
    ###Option=2 -> fastqc of trimmed fastq
    qsub -l h=$Nodo -hold_jid trim_${FastqID} -N QC_trim_${FastqID} -e $Code_Dir/jobs/fastqc.stderr  -o $Code_Dir/jobs/fastqc.stdout $Code_Dir/01_fastqQC.sh $FastqID 2 $Work_Dir/02_trimmed_reads/trimQC $Work_Dir/02_trimmed_reads


    ###############
    # Mapping BWA #
    ###############
    if [ ! -d $Work_Dir/03_mapping_bam ]; then
      mkdir $Work_Dir/03_mapping_bam
      chmod 777 $Work_Dir/03_mapping_bam
    fi
    #Parameters(5): FastqID OutputDir InputDir Option(1) Reference_genome
    ###Option=1 -> fastqc of PE
    ###Option=2 -> fastqc of SE
    qsub -l h=$Nodo -hold_jid QC_trim_${FastqID} -N mapping_${FastqID} -e $Code_Dir/jobs/bwamapping.stderr  -o $Code_Dir/jobs/bwamapping.stdout $Code_Dir/03_MappingBWA.sh $FastqID $Work_Dir/03_mapping_bam $Work_Dir/02_trimmed_reads 1 $Ref_genome


    ##############################
    # Quality control of Mapping #
    ##############################
    if [ ! -d $Work_Dir/03_mapping_bam/mappingQC ]; then
      mkdir $Work_Dir/03_mapping_bam/mappingQC
      chmod 777 $Work_Dir/03_mapping_bam/mappingQC
    fi
    #Parameters(3): FastqID OutputDir InputDir
    qsub -l h=$Nodo -hold_jid $aux -N mapQC_${FastqID} -e $Code_Dir/jobs/mappingQC.stderr  -o $Code_Dir/jobs/mappingQC.stdout $Code_Dir/03_MappingQC.sh $FastqID $Work_Dir/03_mapping_bam/mappingQC $Work_Dir/03_mapping_bam

   
    ######################################################################################################
    # Filtering BAM: filtering unmap, low quality, duplicates, random or unknown chromosomes & blacklist #
    ######################################################################################################
    if [ ! -d $Work_Dir/04_filtering_bam ]; then
      mkdir $Work_Dir/04_filtering_bam
      chmod 777 $Work_Dir/04_filtering_bam
    fi
    #Parameters(4): FastqID OutputDir InputDir BlacklistDir
    qsub -l h=$Nodo -hold_jid mapQC_${FastqID} -N filtering_${FastqID} -e $Code_Dir/jobs/filteringBam.stderr  -o $Code_Dir/jobs/filteringBam.stdout $Code_Dir/04_filteringBam.sh $FastqID $Work_Dir/04_filtering_bam $Work_Dir/03_mapping_bam/mappingQC $BlackList


        ############################
        # Merge BAM (if necessary) #
        ############################
        #Run 04_mergeBam.sh
        #Execution [bash 04_mergeBam.sh SampleNames_Dir Code_Dir Bam_Dir Nodo]

	# MultiQC
	# when all the Fastqc were done, run the multiQC
	# multiqc $Work_Dir/04_merged_filtered_bam -o $Work_Dir/04_merged_filtered_bam/multiQC/multiqc_data/

 
    ###########################
    # Peak Calling Merged BAM #
    ###########################
    if [ ! -d $Work_Dir/05_peak_calling_merged ]; then
      mkdir $Work_Dir/05_peak_calling_merged
      chmod 777 $Work_Dir/05_peak_calling_merged
    fi
    #Parameters(4): FastqID OutputDir InputDir Option(1)
    ###Option=1 -> fastqc of PE (BAMPE)
    ###Option=2 -> fastqc of SE (BAM)
    qsub -l h=$Nodo -hold_jid mergeBam_${FastqID} -N PeakCalling_${FastqID} -e $Code_Dir/jobs/peakCalling.stderr  -o $Code_Dir/jobs/peakCalling.stdout $Code_Dir/05_macs2.sh $FastqID $Work_Dir/05_peak_calling_merged $Work_Dir/04_merged_filtered_bam 1

    
    ######################################
    # Quality Control Peak Calling: FRiP #
    ######################################
    if [ ! -d $Work_Dir/06_frip ]; then
      mkdir $Work_Dir/06_frip
      chmod 777 $Work_Dir/06_frip
    fi  
	#FRiP: Fraction of Reads in Peaks
    #Parameters(4): FastqID InputDir1 InputDir2 OutputDir
    qsub -l h=$Nodo -N frip_${FastqID} -e $Code_Dir/jobs/frip.stderr  -o $Code_Dir/jobs/frip.stdout $Code_Dir/06_FRiP_calculation.sh ${FastqID} $Work_Dir/04_merged_filtered_bam $Work_Dir/05_peak_calling_merged $Work_Dir/06_frip
  

    aux='frip_'${FastqID}
done

 
##############
# FRiP plots #
##############
#Parameters(3): Code_Dir Input_File Work_Dir
qsub -l h=$Nodo -N frip_plot -e $Code_Dir/jobs/frip.stderr  -o $Code_Dir/jobs/frip.stdout $Code_Dir/06_FRiP_plot.sh $Code_Dir frip_input.txt $Work_Dir/06_frip


##########################################################
# Quality Control Peak Calling : Exploration and Jaccard #
##########################################################
if [ ! -d $Work_Dir/07_peak_exploration ]; then
	mkdir /$Work_Dir/07_peak_exploration
	chmod 777 $Work_Dir/07_peak_exploration
fi   
# Jaccard score - Measure of similarity
#Parameters(4): In_File Out_Dir Code_Dir Out_File
qsub -l h=$Nodo -N jaccard -e $Code_Dir/jobs/peakQC.stderr -o $Code_Dir/jobs/peakQC.stdout $Code_Dir/07_peaksQC.sh $Meta_File $Work_Dir/07_peak_exploration $Code_Dir "pairwise_jaccard.txt"


###################
# Consensus Peaks #
###################
if [ ! -d $Work_Dir/08_consensus_peaks ]; then
	mkdir $Work_Dir/08_consensus_peaks
	chmod 777 $Work_Dir/08_consensus_peaks
fi   
#Consensus peaks are generated for each B cell stage and for each individual using GenomicRanges and rtracklayer
#Parameters(3): Code_Dir In_File Work_Dir 
qsub -l h=$Nodo -N consensus -e $Code_Dir/jobs/consensus.stderr  -o $Code_Dir/jobs/consensus.stdout $Code_Dir/08_consensus_peaks.sh $Code_Dir $Meta_File $Work_Dir


##########################################################################
# Quality Control Consensus Peaks: Jaccard score - Measure of similarity #
##########################################################################
#Parameters(3): In_Dir Out_Dir Code_Dir 
qsub -l h=$Nodo -N QCjaccard -e $Code_Dir/jobs/consensusQC.stderr  -o $Code_Dir/jobs/consensusQC.stdout /$Code_Dir/09_consensus_peaksQC.sh $Work_Dir/08_consensus_peaks $Work_Dir/08_consensus_peaks $Code_Dir


################
# Count Matrix #
################
if [ ! -d $Work_Dir/09_count_matrix ]; then
	mkdir $Work_Dir/09_count_matrix
	chmod 777 $Work_Dir/09_count_matrix
fi   
#Parameters(3): In_Dir bed_file Out_Dir
qsub -l h=$Nodo -N count_bed -e $Code_Dir/jobs/bedcounts.stderr  -o $Code_Dir/jobs/bedcounts.stdout $Code_Dir/10_bedCount.sh $Work_Dir/04_merged_filtered_bam $Work_Dir/08_consensus_peaks/consensus_cell_types.bed $Work_Dir/09_count_matrix

#Exploration, normalization and filtering of count matrix.
#Parameters(4): In_Dir Out_Dir Code_Dir meta_file
qsub -l h=$Nodo -N count_bed_QC -e $Code_Dir/jobs/bedcountsQC.stderr  -o $Code_Dir/jobs/bedcountsQC.stdout $Code_Dir/10_bedCount_exploration.sh $Work_Dir/04_merged_filtered_bam $Work_Dir/09_count_matrix $Code_Dir $Meta_File


###################
# Peak Annotation #
###################
if [ ! -d $Work_Dir/10_annotation ]; then
	mkdir $Work_Dir/10_annotation
	chmod 777 $Work_Dir/10_annotation
fi   
#Parameters(3): In_Dir bed_file Out_Dir
qsub -l h=$Nodo -N annotation -e $Code_Dir/jobs/annotation.stderr  -o $Code_Dir/jobs/annotation.stdout $Code_Dir/11_peak_annotation.sh $Work_Dir/08_consensus_peaks/consensus_cell_types.bed $Work_Dir/10_annotation $Work_Dir/09_count_matrix/count_table_with_colnames.txt $Code_Dir


###########################################################
# Differential Binding by transitions + Binary Clustering #
###########################################################
if [ ! -d $Work_Dir/11_diff_binding_transitions ]; then
	mkdir $Work_Dir/11_diff_binding_transitions
	chmod 777 $Work_Dir/11_diff_binding_transitions
fi   
# Differential analysis from normalized and filtered count (TMM) data by cell-stage transitions.
#Parameters(4): Work_Dir counts_file Code_Dir meta_file
qsub -l h=$Nodo -N diff_bind_transitions -e $Code_Dir/jobs/diff_bind_transitions.stderr  -o $Code_Dir/jobs/diff_bind_transitions.stdout $Code_Dir/12_diff_binding_transitions_edgR.sh $Work_Dir/11_diff_binding_transitions $Work_Dir/09_count_matrix/count_table_with_colnames.txt $Code_Dir $Meta_File

#Clustering of significant peaks (TMM data scaled)
if [ ! -d $Work_Dir/11_diff_binding_transitions/clustering ]; then
	mkdir $Work_Dir/11_diff_binding_transitions/clustering
	chmod 777 $Work_Dir/11_diff_binding_transitions/clustering
fi
#Parameters(4): Work_Dir counts_file Code_Dir meta_file
qsub -l h=$Nodo -N transition_peak_clustering -e $Code_Dir/jobs/transition_peak_clustering.stderr  -o $Code_Dir/jobs/transition_peak_clustering.stdout $Code_Dir/13_diff_binding_transitions_clustering.sh $Work_Dir/11_diff_binding_transitions/sig_dea_transitions.txt $Work_Dir/10_annotation/consensus_peaks_annotation.txt $Work_Dir $Code_Dir    
