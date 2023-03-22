#!/bin/bash
#SBATCH --job-name=peak_calling
#SBATCH --output=peak_calling.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=4
#SBATCH --mem=90G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl


# this script is to call peak per group with the approriate control and everything together without control 
# Nhung, 21 03 2023

# Tool: 
# - macs2 ver 2.2.7.1


# Input file options

# -t: The IP data file (this is the only REQUIRED parameter for MACS)
# -c: The control or mock data file
# -f: format of input file; Default is “AUTO” which will allow MACS to decide the format automatically.
# -g: mappable genome size which is defined as the genome size which can be sequenced; some precompiled values provided.
# Output arguments
# --outdir: MACS2 will save all output files into speficied folder for this option
# -n: The prefix string for output files

# Input: 
# - per group: 1 bam file for the test, 1 bam file for the control (after removing duplicate and filtering)
# 
# Output:
# -  broad: .broadPeak (most important), .gappedPeak, .xls 
# -  narrow: .narrowPeak (most important), .xls, .bed
######### Define global variables ################

data_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/SCC_ChIC-PMC-DRO_plates_20210520_run1

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

peak_dir=${res_dir}/peakCalling/${sample_ID}
mkdir -p ${peak_dir}

peak_no_control_dir=${res_dir}/peakCalling_nocontrol
mkdir -p ${peak_no_control_dir}

# tool dir


bowtie2Index=/hpc/pmc_drost/SOURCES/Genomes/human/bowtie2/human_gencode37_hg38

picardTool=/hpc/pmc_drost/PROJECTS/swang/software/picard.jar
new_tmp_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/tmp # to solve the out of space with the temporary output from picard



# sample_Ids was generated as text file from ls filename.txt from the data_dir

sample_IDs=( "bulkChIC-PMC-DRO-011" \
            "bulkChIC-PMC-DRO-012"\
            "bulkChIC-PMC-DRO-013"\
            "bulkChIC-PMC-DRO-014"\
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
            "SCC-ChIC-PMC-DRO-L1")

total_sample=${#sample_IDs[@]}

# classify sample for peakcalling
tfe3=("SCC-ChIC-PMC-DRO-T1 SCC-ChIC-PMC-DRO-T5 bulkChIC-PMC-DRO-016 SCC-bulkChIC-PMC-DRO-008")
luciferase=("SCC-ChIC-PMC-DRO-L1 SCC-ChIC-PMC-DRO-L5 bulkChIC-PMC-DRO-014 SCC-bulkChIC-PMC-DRO-005")
fusion=("SCC-ChIC-PMC-DRO-F1 SCC-ChIC-PMC-DRO-F5 bulkChIC-PMC-DRO-015 SCC-bulkChIC-PMC-DRO-002")
allT=("$tfe3 $luciferase $fusion")
tfe3C=$res_dir/rm_dup/bulkChIC-PMC-DRO-013/bulkChIC-PMC-DRO-013_rmdup_filt.bam
luciferaseC=$res_dir/rm_dup/bulkChIC-PMC-DRO-011/bulkChIC-PMC-DRO-011_rmdup_filt.bam
fusionC=$res_dir/rm_dup/bulkChIC-PMC-DRO-012/bulkChIC-PMC-DRO-012_rmdup_filt.bam


# Run per group using appropriate controls 
task () {
  output_dir=${peak_dir}/${sample_ID}
  mkdir -p ${output_dir}

  sample_dir=( $(find ${rm_dup_dir}/${sample_ID} -name "*rmdup_filt.bam") )
  
  echo "----------------start-narrow-peak-calling-for $sample_ID at $(date)----------"
  macs2 callpeak --tempdir $new_tmp_dir -t ${sample_dir} -c ${control} -f BAMPE -g hs -n ${sample_ID}_paired_control -q 0.01 --outdir $output_dir/narrow

  echo "----------------start-broad-peak-calling-for $sample_ID at $(date)----------"

  macs2 callpeak --tempdir $new_tmp_dir -t ${sample_dir} -c ${control} -f BAMPE -g hs -n ${sample_ID}_paired_control -q 0.01 --broad --outdir $output_dir/broad

  echo "----------------finish-broad-peak-calling-for $sample_ID at $(date)----------" ;

}

# run for tfe3
total_sample=${#tfe3[@]}
n=0

for sample_ID in $tfe3; do
  n=$((n+1))
  echo "-----------running $n out of $total_sample tfe3 samples---------------"
  
  control=$tfe3C

  task "$sample_ID" &
  
done
wait
echo "all done for tfe3"

# run for luciferase
total_sample=${#luciferase[@]}
n=0
for sample_ID in $luciferase;do
  n=$((n+1))
  echo "-----------running $n out of $total_sample luciferase samples---------------"

  control=$luciferaseC

  task "$sample_ID" &

done

wait

echo "all done for luciferase"

# run for fusion

total_sample=${#fusion[@]}
n=0
for sample_ID in $fusion;do
  n=$((n+1))
  echo "-----------running $n out of $total_sample fusion samples---------------"

  control=$fusionC

  task "$sample_ID" &

done

wait

echo "all done for fusion"

# Run for all samples without control

task () {
  output_dir=${peak_no_control_dir}/${sample_ID}
  mkdir -p ${output_dir}

  sample_dir=( $(find ${rm_dup_dir}/${sample_ID} -name "*rmdup_filt.bam") )
  
  echo "----------------start-narrow-peak-calling-no-control-for $sample_ID at $(date)----------"
  macs2 callpeak --tempdir $new_tmp_dir -t ${sample_dir} -f BAMPE -g hs -n ${sample_ID}_paired -q 0.01 --outdir $output_dir/narrow

  echo "----------------start-broad-peak-calling-for $sample_ID at $(date)----------"

  macs2 callpeak --tempdir $new_tmp_dir -t ${sample_dir} -f BAMPE -g hs -n ${sample_ID}_paired -q 0.01 --broad --outdir $output_dir/broad

  echo "----------------finish-broad-peak-calling-no-control-for $sample_ID at $(date)----------" ;

}


total_sample=${#allT[@]}
n=0
for sample_ID in $allT;do
  n=$((n+1))
  echo "-----------running $n out of $total_sample all samples---------------"

  control=$fusionC

  task "$sample_ID" &

done

wait

echo "all done for no control"


