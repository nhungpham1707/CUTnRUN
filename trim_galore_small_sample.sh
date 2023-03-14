#!/bin/bash
#SBATCH --job-name=trim_galore_small_test
#SBATCH --output=test_small_samples_galore.out
#SBATCH --time=10:0:0
#SBATCH -N 3
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
# data dir
data_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/SCC_ChIC-PMC-DRO_plates_20210520_run1


res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test


# sample_Ids was generated as text file from ls filename.txt from the data_dir
sample_IDs=( "bulkChIC-PMC-DRO-011" \
			"bulkChIC-PMC-DRO-013")


# find read IDs from a sample

for str in ${sample_IDs[@]}; do
  
  sample_ID=$str
  echo $sample_ID
    
  fastq_IDs=( $(find ${data_dir}/$sample_ID -name "*.fastq.gz") )
  #declare -p fastq_IDs

  len=${#fastq_IDs[@]}
  echo $len
  echo "start trimming for $sample_ID at $(date)"  	

 # for (( i=0; i<$len; i++ )); do echo "${fastq_IDs[$i]}" ; done
# #Adapter and quality trimming (Trim Galore!)
  	mkdir -p ${res_dir}/trimmed/${sample_ID}
 	sample_trim=${res_dir}/trimmed/${sample_ID}
 	trim_galore --fastqc --paired ${fastq_IDs[0]} ${fastq_IDs[1]} -o $sample_trim
	
	echo "finish trimming for $sample_ID at $(date)"
	# fastqc SCC-bulkChIC-PMC-DRO-005_AAALH7CM5_S8_L001${suffixF} -o $basedir
	# fastqc SCC-bulkChIC-PMC-DRO-005_AAALH7CM5_S8_L001${suffixR} -o $basedir
done

