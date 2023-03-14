#!/bin/bash
#SBATCH --job-name=fastqc_all_samples
#SBATCH --output=before_trimming_fastqc_all_samples.out
#SBATCH --time=96:0:0
#SBATCH -N 3
#SBATCH --mem=90G
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl


# this script is to conduct fastqc prior to trimming for all samples 

# Nhung, 13 03 2023

# Tool: fastqc v0.12.1

# Input: R1 and R2 fastq sequences for each sample

# Output:
# - output from each sample will be saved in a subfolder with the respective sample ID
# - For each sample a zip file and a html file

#########
# data dir
data_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/SCC_ChIC-PMC-DRO_plates_20210520_run1


res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test

# sample_Ids was generated as text file from ls filename.txt from the data_dir
sample_IDs=( "bulkChIC-PMC-DRO-011" \
      "bulkChIC-PMC-DRO-013"\
      "bulkChIC-PMC-DRO-016"\
      "SCC-bulkChIC-PMC-DRO-005"\
      "SCC-ChIC-PMC-DRO-L5"\
      "bulkChIC-PMC-DRO-012"\
      "SCC-bulkChIC-PMC-DRO-008"\
      "SCC-ChIC-PMC-DRO-LH"\
      "bulkChIC-PMC-DRO-013"\
      "SCC-ChIC-PMC-DRO-F1"\
      "SCC-ChIC-PMC-DRO-T1"\
      "SCC-ChIC-PMC-DRO-F5"\
      "SCC-ChIC-PMC-DRO-T5"\
      "bulkChIC-PMC-DRO-014"\
      "bulkChIC-PMC-DRO-015"\
      "SCC-ChIC-PMC-DRO-FH"\
      "SCC-ChIC-PMC-DRO-TH"\
      "SCC-bulkChIC-PMC-DRO-002"\
      "SCC-ChIC-PMC-DRO-L1")


# find read IDs from a sample

for str in ${sample_IDs[@]}; do
  
  sample_ID=$str
  echo $sample_ID
    
  fastq_IDs=( $(find ${data_dir}/$sample_ID -name "*.fastq.gz") )
  #declare -p fastq_IDs

  len=${#fastq_IDs[@]}
  echo $len
 # for (( i=0; i<$len; i++ )); do echo "${fastq_IDs[$i]}" ; done
# #Adapter and quality trimming (Trim Galore!)
 	#cd $trimmed
 	mkdir -p ${res_dir}/fastqc_before_trimming/${sample_ID}
 	sample_fastqc=${res_dir}/fastqc_before_trimming/${sample_ID}
  echo "start fastqc for $sample_ID at $(date)"

  fastqc ${fastq_IDs[0]} -o $sample_fastqc
  fastqc ${fastq_IDs[1]} -o $sample_fastqc

  echo "finish fastqc for $sample_ID at $(date)"
done


