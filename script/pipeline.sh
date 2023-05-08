#!/bin/bash
#SBATCH --job-name=diffbind
#SBATCH --output=diffbind.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=90G
#SBATCH --cpus-per-task=8
#SBATCH --gres=tmpspace:30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl

#==================================================================
#   This is a master script to analyze cut and run data. 
#   The script can also be adapted for ATAC, ChiPSeq or 
#   data that required similar steps. 
#
#   More infomation can be found at 
#   https://github.com/nhungpham1707/CUTnRUN
#   Author: Nhung, 01 05 2023
#==================================================================

#==================================================================
#                   TOOLS & FLAGS EXPLAINATION
#
#   Quality check:
#   - fastqc v0.12.1
#
#   Trimming: 
#   - trim_galore ver 0.6.10

#   Alignment: 
#	- bowtie ver 2.5.1, input flags:
#			-end-to-end: for trimmed sequence
#       (used local if seq is not trimmed before this step)
#			-p 16: 16 threads
#			-samtools -bS: to save input sam file as bam file
# 
#   Remove dup and filtering for short reads:
#   - picard version 52.0  
#   (java -jar /hpc/pmc_drost/PROJECTS/swang/software/picard.jar   MarkDuplicates —version)
#   - samtools, ver 1.16.1, input flags:
# 	(https://wikis.utexas.edu/display/CoreNGSTools/Filtering+with+)
# 	(http://www.htslib.org/doc/samtools-view.html)
# 		-b to get output in bam file format
# 		-q 20 skip alignment smaller than 20
# 		-F do not output alignment with any bit set in flag 
#       (1024 mean from pcr or optical dups)
# 		-h keep header
# 		-f only output alignment with all bits set in flag 
#    (2 = read map in proper pair https://www.samformat.info/sam-format-flag)
#
#   Convert bam file to bigwig and bedgraph
#   - bamcoverage 3.5.1   (/hpc/pmc_drost/nhung/anaconda3/envs/cutnrun_trimgalore/bin/bamCoverage --version)
#   - Input flags:
#     -b input files
#     -effectiveGenomeSize  2913022398 for hg38. 
#       for other genomes can be found here ref. https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html
#     -ignoreForNormalization  indicate chromosome that will not be used for normalization
#     -o output 
#     - extra flag for bedgraph
#     --outFileFormat bedgraph
#     --binSize $bin_size # uncomment out lines: L250-L258 writeBedGraph.py 
#       in deeptools if want to generate files with a fix bin size
#
#   Data normalization:
#   - s3norm v1
#   follow set up here https://github.com/guanjue/S3norm.git 
#   or create the conda env with spec-file-s3norm.txt and clone the github
#
#   Peak calling
#   - macs2 ver 2.2.7.1 
#   - input flags:
#        -t: The IP data file (this is the only REQUIRED parameter for MACS)
#        -c: The control or mock data file 
#       (if not specify macs2 will use a general genome as control)
#        -f: format of input file; 
#   Default is “AUTO” which will allow MACS to decide the format automatically. In BAMPE mode, there is no need to estimate fragment lengths since the actual insertion length for each read pair will be considered. Therefore no model.R will be generated. https://github.com/macs3-project/MACS/issues/281
#        -g: mappable genome size which is defined as the genome size which can be sequenced; some precompiled values provided.
#       Output arguments
#       --outdir: MACS2 will save all output files into speficied folder for this option
#       -n: The prefix string for output files
#       -q: The false discovery rate cut-off. default 0.01 


# diffBind version 2.28
#==================================================================

#==================================================================
#                       INPUT & OUTPUT DATA
#
#   Quality check 
#       Input: R1 and R2 fastq sequences for each sample
#       Output:
#       - output from each sample will be saved in a subfolder with the respective sample ID
#       - For each sample a zip file and a html file
#
#   Trimming 
#       Input: R1 and R2 fastq sequences for each sample
#       Output:
#           - 2 trimming reports for R1 and R2.
#           - 2 val_ files for R1 and R2 (inputs for next step)
#           - 2 fastqc html for R1 and R2
#           - 2 fastqc.zip files for R1 and R2.
#
#   Alignment
#       Input: R1 and R2 outputs from trimming for each sample
#       Output:
#           - 1 bam file per sample stored in its respective folder 
#
#   Remove dup and filtering for short reads: 
#		Input data: 1 bam file for each sample (output from alignment)
# 	    Output:
# 		- 1 mark and remove duplicate bam file
# 		- 1 filter bam file (input for next step peakcalling)
# 		- 1 index bam.bai file
# 		- 1 metric text file
#
#   Data normalization with s3norm 
#       Input data: 
#       - bedgraph files with fix bin size and no 0 values. 
#       (The input files can be generated from bam file after removing dups in modifybedgraph.sh)
#       - a sample_info.txt file. Make sure the columns are proper tab delimited and keep an empty line at the end of the file
#       Outputs: 
#       - S3norm_rc_bedgraph: normalize counts. Use for EdgeR, DESeq2 or differential peak calling
#       - NBP_bedgraph: The negative log10 p-value of S3norm normalized read counts based on a negative binomial background model
#       - S3norm_NBP_bedgraph: The S3norm normalized negative log10 p-value based on a negative binomial background model. These files are used for peakcalling 'bdgpeakcall' and 'bdgbroadcall' in MACS2 (https://github.com/taoliu/MACS))
#
#   Peak calling before data normalization 
#       Inputs: 
#       With control sample
#       - per group: 1 bam file for the test, 
#       1 bam file for the control (after removing duplicate and filtering)
#       Without control sample
#        - all bam files after removing dup and filtering
#
#       Outputs:
#       -  broad: .broadPeak (most important), .gappedPeak, .xls 
#       -  narrow: .narrowPeak (most important), .xls (input for diffBind), summits.bed (input for heatmap) 
#       summits.bed has information of a peak where the score (protein binding intensity) is maximum in the peak whereas the narrowPeak file has complete peak boundary. The summits.bed reports a peak with width/distance 1bp. That means this is the region (of peak) where the protein binding intensity reach its maximum. Therefore every peak should have its own summit. ref. https://www.biostars.org/p/9470414/
# more about summit and extraction of seq can be found here https://notebook.community/ssjunnebo/pathogen-informatics-training/Notebooks/ChIP-Seq/motif-analysis

#   Peakcalling after data normalization
#       Inputs:
#       - bedgraph files in S3norm_NBP_bedgraph (change file extension name in 6b-peakCalling_normalize.sh if needed)
#       Outputs:
#       - narrow peak files
#==================================================================

#==================================================================
#                   DEFINE GLOBAL VARIABLES
#
# source /hpc/pmc_drost/nhung/anaconda3/ect/profile.d/conda.sh
# conda activate cutnrun_trimgalore
# Adap data_dir and res_dir before running
data_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/SCC_ChIC-PMC-DRO_plates_20210520_run1
# dir for new data 
# data_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/SCC_ChIC-PMC-DRO_plates_20210520_run1/SCC-bulkChIC-PMC-DRO-020--25/
# result dir
res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test

qc_dir=${res_dir}/qualityCheck
mkdir -p $qc_dir

trim_dir=${res_dir}/trimmed
mkdir -p $trim_dir

align_dir=${res_dir}/alignment
mkdir -p $align_dir

rm_dup_dir=$res_dir/rm_dup
mkdir -p $rm_dup_dir

frip_dir=${res_dir}/FRiP
mkdir -p $frip_dir

s3norm_script_directory='/hpc/pmc_drost/nhung/S3norm'
s3norm_working_directory=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/modify_bedgraph/

peak_dir=${res_dir}/peakCalling
mkdir -p ${peak_dir}

peak_no_control_dir=${res_dir}/peakCalling_nocontrol
mkdir -p ${peak_no_control_dir}

normalize_peak_dir=${res_dir}/normalize_peakCalling
mkdir -p ${normalize_peak_dir}

motif_dir=${res_dir}/motif
mkdir -p ${motif_dir}

merged_bigwig_dir=${res_dir}/merged_bigwig
mkdir -p ${merged_bigwig_dir}

bigwig_dir=${res_dir}/bigwig
mkdir -p $bigwig_dir

peak_analysis_dir=${res_dir}/peak_analysis
mkdir -p $peak_analysis_dir

figure_dir=${res_dir}/figures
mkdir -p $figure_dir

bincount_dir=${res_dir}/bincount_window
mkdir -p $bincount_dir

diffBind_res_dir=${res_dir}/diffBind_analysis
mkdir -p $diffBind_res_dir
# tool dir

bowtie2Index=/hpc/pmc_drost/SOURCES/Genomes/human/bowtie2/human_gencode37_hg38

picardTool=/hpc/pmc_drost/PROJECTS/swang/software/picard.jar
new_tmp_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/tmp # to solve the out of space with the temporary output from picard

homer_dir=/hpc/pmc_drost/nhung/anaconda3/envs/cutnrun_trimgalore/bin/homer
findMotif_dir=/hpc/pmc_drost/nhung/anaconda3/envs/cutnrun_trimgalore/bin/findMotifsGenome.pl

bamCoverage_dir=/hpc/pmc_drost/nhung/anaconda3/envs/cutnrun_trimgalore/bin/bamCoverage # to merge 
effectiveGenomeSize=2913022398 
# fasta genome directory for motif finding. It has to be the same with the reference genome that was used for alignment
fasta_genome_dir=/hpc/pmc_drost/SOURCES/Genomes/human/gencode37_GRCh38_primary_assembly_genome.fa

# hg38 genes list directory to generate heatmap
hg38_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/hg38_gene_2.bed 
 
# sample_Ids was generated as text file from ls filename.txt from the data_dir

# sample_IDs=( "bulkChIC-PMC-DRO-011" \
#             "bulkChIC-PMC-DRO-012"\
#             "bulkChIC-PMC-DRO-013"\
#             "bulkChIC-PMC-DRO-014"\
#             "bulkChIC-PMC-DRO-015"\
#             "bulkChIC-PMC-DRO-016"\
#             "SCC-bulkChIC-PMC-DRO-002"\
#             "SCC-bulkChIC-PMC-DRO-005"\
#             "SCC-bulkChIC-PMC-DRO-008"\
#             "SCC-ChIC-PMC-DRO-L5"\
#             "SCC-ChIC-PMC-DRO-LH"\
#             "SCC-ChIC-PMC-DRO-F1"\
#             "SCC-ChIC-PMC-DRO-F5"\
#             "SCC-ChIC-PMC-DRO-FH"\
#             "SCC-ChIC-PMC-DRO-T1"\
#             "SCC-ChIC-PMC-DRO-T5"\
#             "SCC-ChIC-PMC-DRO-TH"\ 
#             "SCC-ChIC-PMC-DRO-L1"\ 
#             "SCC-bulkChIC-PMC-DRO-020"\
#             "SCC-bulkChIC-PMC-DRO-021"\
#             "SCC-bulkChIC-PMC-DRO-022"\
#             "SCC-bulkChIC-PMC-DRO-023"\
#             "SCC-bulkChIC-PMC-DRO-024"\
#             "SCC-bulkChIC-PMC-DRO-025")
# sample_IDs=( "bulkChIC-PMC-DRO-011" "SCC-ChIC-PMC-DRO-T1")
#
# sample IDs for peakCalling after normalization (remove samples that are used as control for data normalization)
sample_IDs=( "bulkChIC-PMC-DRO-014"\
            "bulkChIC-PMC-DRO-015"\
            "bulkChIC-PMC-DRO-016"\
            "SCC-bulkChIC-PMC-DRO-002"\
            "SCC-bulkChIC-PMC-DRO-005"\
            "SCC-bulkChIC-PMC-DRO-008"\
            "SCC-ChIC-PMC-DRO-L5"\
            "SCC-ChIC-PMC-DRO-LH"\
            "SCC-ChIC-PMC-DRO-F1"\
            "SCC-ChIC-PMC-DRO-F5"\
            "SCC-ChIC-PMC-DRO-FH"\
            "SCC-ChIC-PMC-DRO-T1"\
            "SCC-ChIC-PMC-DRO-T5"\
            "SCC-ChIC-PMC-DRO-TH"\ 
            "SCC-ChIC-PMC-DRO-L1"\ 
            "SCC-bulkChIC-PMC-DRO-020"\
            "SCC-bulkChIC-PMC-DRO-021"\
            "SCC-bulkChIC-PMC-DRO-022"\
            "SCC-bulkChIC-PMC-DRO-023"\
            "SCC-bulkChIC-PMC-DRO-024"\
            "SCC-bulkChIC-PMC-DRO-025")
            
total_sample=${#sample_IDs[@]}

# classify sample for peakcalling
tfe3=( "SCC-ChIC-PMC-DRO-T1" "SCC-ChIC-PMC-DRO-T5" "bulkChIC-PMC-DRO-016" "SCC-bulkChIC-PMC-DRO-008")
luciferase=( "SCC-ChIC-PMC-DRO-L1" "SCC-ChIC-PMC-DRO-L5" "bulkChIC-PMC-DRO-014" "SCC-bulkChIC-PMC-DRO-005")
fusion=( "SCC-ChIC-PMC-DRO-F1" "SCC-ChIC-PMC-DRO-F5" "bulkChIC-PMC-DRO-015" "SCC-bulkChIC-PMC-DRO-002")
# new luc and fusion samples
# luciferase="SCC-bulkChIC-PMC-DRO-020 SCC-bulkChIC-PMC-DRO-021 SCC-bulkChIC-PMC-DRO-022"
# fusion="SCC-bulkChIC-PMC-DRO-023 SCC-bulkChIC-PMC-DRO-024 SCC-bulkChIC-PMC-DRO-025"

#allT="$tfe3 $luciferase $fusion"
tfe3C=$res_dir/rm_dup/bulkChIC-PMC-DRO-013/bulkChIC-PMC-DRO-013_rmdup_filt.bam
luciferaseC=$res_dir/rm_dup/bulkChIC-PMC-DRO-011/bulkChIC-PMC-DRO-011_rmdup_filt.bam
fusionC=$res_dir/rm_dup/bulkChIC-PMC-DRO-012/bulkChIC-PMC-DRO-012_rmdup_filt.bam
h3k4me3=( "SCC-ChIC-PMC-DRO-FH" "SCC-ChIC-PMC-DRO-LH" "SCC-ChIC-PMC-DRO-TH" )
h3k4me3_2=( "bulkChIC-PMC-DRO-011" "bulkChIC-PMC-DRO-012" "bulkChIC-PMC-DRO-013") 
h3k4me3_all=( "SCC-ChIC-PMC-DRO-FH" "SCC-ChIC-PMC-DRO-LH" "SCC-ChIC-PMC-DRO-TH" "bulkChIC-PMC-DRO-011" "bulkChIC-PMC-DRO-012" "bulkChIC-PMC-DRO-013")

# define what samples to find motif 
motif_samples=$allT
#==================================================================

#==================================================================
#                           RUNNING ANALYSIS
# turn of core file generation to not overload home dir
 ulimit -c 0 # to disable coredump  https://community.hpe.com/t5/system-administration/how-to-disable-or-restrict-core-dumps/td-p/4276097#.ZDzoIpNByDU

echo "start running cut and run analysis at $(date)"
#==================================================================

#==================================================================
#                          BEFORE PEAK CALLING
# step 1. sequencing quality check 
# echo "------------------step1. running quality check--------------"
# . ./1-qualityCheck.sh

# # # step 2. adapter and bad reads trimming 
# echo "-------------------step 2. running trimming-----------------"
# . ./2-trimming.sh 

# # # step 3. alignment- map to hg38 genome 
# echo "-------------------step 3. running alignment----------------"
# . ./3-alignment.sh

# # # step 4. filtering: remove duplciates and reads < 20bp
# echo "-------------------step 4. running filtering----------------"
# . ./4-filtering.sh

# # step 5. merge and transform bam file to bigwig
# echo "-------------------step 5. running transform bam to bigwig--"
# . ./5-bam2bigwig.sh

# step 6. Calculate fraction of read in peak (FRiP)
# echo "-------------------step 6. running frip calculation-------- "
# . ./6-Calculate_FRiP.sh

# step 7. Normalize data
# echo "-------------------step 7. running data normalization------ "
# s3norm_sample_file_name='no_zeros_testsamples.txt'
#. ./7-run_s3norm.sh 
#==================================================================

#==================================================================
#                   PEAK CALLING
# # step 8a. peak calling with raw data
# echo "---------step 8a. running peak calling for raw data---------"
# . ./8a-peakCalling.sh

# step 8b. peak calling with normalize data
# echo "-------------------step 8. running peak calling after normalization----------------------"
# . ./8b-peakCalling_normalize.sh
#==================================================================

#==================================================================
#                   AFTER PEAKCALLING

# step 9. Extract overlap peak from replicates in the same condition. # modify variable names if using for different experiments before run
# echo "-------------------step 9. running peak overlap---------------"
# . ./9-GetPeakOverlap.sh # few minutes 

# step 10. Extract peak overlap statistic
#  echo "-------------------step 10. running peak statistic ---------------"
#  Rscript 10-Peak_statistic.R

# step 11. Identify peaks that are differentially enriched between conditions. Modify variable names if using for different experiments before run
echo "---------------step 11. running peak differential analysis---------------"
export DIFFBIND_RESULT_DIR_VARIABLE=$diffBind_res_dir
export SAMPLE_SHEET_DIR_VARIABLE=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/R/diffbind_normalize_samples.csv
Rscript DiffBind_analysis.R

# step 12. heatmap generation. prior to run: change sample paths in 9-heatmap.sh to those that one wish to make the heatmap for and if require also the bed file that indicate the desire genome region to plot. 
# echo "---------------step 12. running heatmap generation---------------"
# . ./12-heatmap.sh
. ./12-heatmap.sh "fusionvshg38" "$hg38_dir" "${fusion[@]}"

# step 13. prepare for motif analysis. Prior to run change sample paths in 10-prepareMotifAnalysis.sh if needed. 
# echo "--------------------step 13. running motif finding preparation"
# . ./13-prepareMotifAnalysis.sh

# step 14. motif finding 
# echo "-------------------step 14. running motif finding----------------------"
# . ./14-motifFinding.sh 

# step 15. motif annotation 
# step 16. peak annnotation
# step 17. super enhancer finding 
# step. check sample correlation
# . ./SamplesCorrelation.sh

# convert bam file to bed
# . ./bam2bed.sh
# . ./binCount.sh


##### Data normalization ###############
## prepare bedgraph files with the same bin size for all samples
# . ./bincount.sh # need to test

## modify bedgraph file to remove rows with 0 count in all samples. If 0 values are more than 10% in the sample, s3norm will fail (log0 is inf), hence add 1 to these 0 values in all samples
# python modifyBedgraphForS3norm.py # already test - worked!

# start normalization on the modified bedgraph files



# . ./run_s3norm.sh
# calculate frip 
# . ./Calculate_FRiP.sh

#==================================================================

#==================================================================
#                           GENERATE PLOTS AND REPORTS

# step 7. Generate overview plots to access number of reads, duplicate, frip and peaks 
# echo "-------------------step 7. running report plots----------------------"
# Rscript 7-Report_plots.R

# step 8. generate plots from histone samples to check the reliability of the cut and run experiment


echo "finish cut and run analysis at $(date)"

