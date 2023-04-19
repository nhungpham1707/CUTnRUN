#!/bin/bash
#SBATCH --job-name=alignment_all
#SBATCH --output=alignment_all.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=90G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl


# this script is to run alignment for all samples
# Nhung, 21 03 2023

# Tool: 
# - bowtie ver 0.6.10 , input flags:
#			-end-to-end: for trimmed sequence
#			-p 16: 16 threads
#			-samtools -bS: to save input sam file as bam file

# Input: R1 and R2 outputs from trimming for each sample

# Output:
# - 1 bam file per sample stored in its respective folder


      
########### Make task function ##############
task () {
trim_IDs=( $(find ${trim_dir}/${sample_ID} -name "*.fq.gz") )

len=${#trim_IDs[@]}
  
echo $len

sample_align_dir=${align_dir}/${sample_ID}
mkdir -p $sample_align_dir

 echo "-----start alignment for $sample_ID at $(date)-------"

# Alignment for trimmed seq 

bowtie2 --end-to-end --very-sensitive --no-mixed --no-discordant --dovetail -p 16 -x $bowtie2Index -1 ${trim_IDs[0]} -2 ${trim_IDs[1]} -S ${sample_align_dir}/${sample_ID}.sam  &> ${sample_align_dir}/${sample_ID}_bowtie2.txt  

echo "-------finish alignment for $sample_ID at $(date)--------" ;
}

#### run loop

n=0

for sample_ID in ${sample_IDs[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample samples---------------------- "

    task '$sample_ID' &

done

wait
echo "all done"


     
