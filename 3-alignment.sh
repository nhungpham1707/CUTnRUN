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

bowtie2 --end-to-end -p 16 -x $bowtie2Index -1 ${trim_IDs[0]} -2 ${trim_IDs[1]} 1>${sample_align_dir}/${sample_ID}.sam 2> ${sample_align_dir}/${sample_ID}.txt 

samtools view -bS ${sample_align_dir}/${sample_ID}.sam > ${sample_align_dir}/${sample_ID}.bam

# extract the fragment length from the 9th column to make report figure (check report.r script)
samtools view -F 0x04 ${sample_align_dir}/${sample_ID}.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' >${sample_align_dir}/${sample_ID}_fragmentLen.txt

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


     
