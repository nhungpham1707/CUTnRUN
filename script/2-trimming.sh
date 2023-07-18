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

mkdir -p $trim_dir

#########
task () {

  fastq_IDs=( $(find ${data_dir}/$sample_ID -name "*.$input_extension") )
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
  trim_file=( $(find ${trim_dir}/${sample_ID} -maxdepth 1 -name "*.fq.gz") )
  if [ -s "$trim_file" ] 
  then
    echo "$trim_file exists. Skip trimming step for $sample_ID!"
  else 
    task '$sample_ID' &
  fi
  done

  wait
multiqc_file=$trim_dir/trim_report.html
 if [ -s "$multiqc_file" ] 
  then
    echo "$multiqc_file exists. Skip generating qc report for trimming step!"
  else 
  # generate report
  multiqc ${trim_dir} -n trim_report -o ${trim_dir}
fi
  echo "all done trimming"
