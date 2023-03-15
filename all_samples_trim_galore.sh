#!/bin/bash
#SBATCH --job-name=trim_galore_all_samples
#SBATCH --output=all_samples_galore.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=90G
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl


# this script is to conduct trimming for all samples 

# Nhung, 09 03 2023

# Tool: trim_galore, ver 0.6.10

# Input: R1 and R2 fastq sequences for each sample

# Output:
# - 2 trimmed files for R1 and R2.
# - 2 trimming reports for R1 and R2.
# - 2 val_ files for R1 and R2. 
# - 2 fastqc html for R1 and R2 

#########

#conda activate cutnrun_trimgalore

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



for str in ${sample_IDs[@]}; do
  
  sample_ID=$str
  echo $sample_ID
    
  fastq_IDs=( $(find ${data_dir}/$sample_ID -name "*.fastq.gz") )
  #declare -p fastq_IDs

  len=${#fastq_IDs[@]}
  echo $len
 # for (( i=0; i<$len; i++ )); do echo "${fastq_IDs[$i]}" ; done
# #Adapter and quality trimming (Trim Galore!)
 	mkdir -p ${res_dir}/trimmed/${sample_ID}
 	sample_trim=${res_dir}/trimmed/${sample_ID}
  
  echo "start trim_galore for $sample_ID at $(date)"

 	trim_galore --fastqc --paired ${fastq_IDs[0]} ${fastq_IDs[1]} -o $sample_trim

  echo "finish trim_galore for $sample_ID at $(date)"


done


