#!/bin/bash
#SBATCH --job-name=paralel_fastqc_all_samples
#SBATCH --output=paralel_fastqc.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=3
#SBATCH --mem=40G
#SBATCH --cpus-per-task=5
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
      "SCC-bulkChIC-PMC-DRO-005"\)


# find read IDs from a sample

total_sample=${#sample_IDs[@]}
n=0
for str in ${sample_IDs[@]}; do
  n=$((n+1))
  echo "-----------running $n out of $total_sample samples---------------------- "
  
  sample_ID=$str
  echo $sample_ID
    
  fastq_IDs=( $(find ${data_dir}/$sample_ID -name "*.fastq.gz") )
  #declare -p fastq_IDs

  len=${#fastq_IDs[@]}
  echo $len
 # for (( i=0; i<$len; i++ )); do echo "${fastq_IDs[$i]}" ; done
  mkdir -p ${res_dir}/fastqc_before_trimming/${sample_ID}
  sample_fastqc=${res_dir}/fastqc_before_trimming/${sample_ID}
  
  echo "start fastqc for $sample_ID at $(date)"

  fastqc ${fastq_IDs[0]} -o $sample_fastqc
  fastqc ${fastq_IDs[1]} -o $sample_fastqc

  echo "finish fastqc for $sample_ID at $(date)" &

done

wait

echo("all done")



