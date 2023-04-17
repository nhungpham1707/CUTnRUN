#!/bin/bash
#SBATCH --job-name=trimming
#SBATCH --output=trimming.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=90G
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl


# this script is to conduct trimming for all samples

# Nhung, 21 03 2023

# Tool: trim_galore, ver 0.6.10

# Input: R1 and R2 fastq sequences for each sample

# Output:
# - 2 trimming reports for R1 and R2.
# - 2 val_ files for R1 and R2 (inputs for next step)
# - 2 fastqc html for R1 and R2
# - 2 fastqc.zip files for R1 and R2.

#########
task () {

  fastq_IDs=( $(find ${data_dir}/$sample_ID -name "*.fastq.gz") )
  #declare -p fastq_IDs

  len=${#fastq_IDs[@]}
  echo $len
  
  sample_trim=${trim_dir}/${sample_ID}
  mkdir -p ${sample_trim}
        
  echo "start trim_galore for $sample_ID at $(date)"

  trim_galore --fastqc --paired ${fastq_IDs[0]} ${fastq_IDs[1]} -o $sample_trim

  echo "finish trim_galore for $sample_ID at $(date)" ;
}

n=0

for sample_ID in ${sample_IDs[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample samples---------------------- "

    task '$sample_ID' &

done

wait

# generate report
multiqc ${trim_dir}
echo "all done trimming"

