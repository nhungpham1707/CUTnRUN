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
# Nhung, 20 03 2023

# Tool: 
# - bowtie ver 0.6.10 , input flags:
#			-end-to-end: for trimmed sequence
#			-p 16: 16 threads
#			-samtools -bS: to save input sam file as bam file

# Input: R1 and R2 outputs from trimming for each sample

# Output:
# - 1 bam file per sample stored in its respective folder

######### define global variables ################

data_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/SCC_ChIC-PMC-DRO_plates_20210520_run1

# result dir
res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test
align_dir=${res_dir}/alignment
mkdir -p $align_dir

# tool
bowtie2Index=/hpc/pmc_drost/SOURCES/Genomes/human/bowtie2/human_gencode37_hg38

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
n=0
      
########### Make task function ##############
task () {
trim_IDs=( $(find ${res_dir}/trimmed/${sample_ID} -name "*.fq.gz") )

len=${#trim_IDs[@]}
  
echo $len

sample_align_dir=${align_dir}/${sample_ID}
mkdir -p $sample_align_dir

 echo "-----start alignment for $sample_ID at $(date)-------"

# Alignment for trimmed seq 

bowtie2 --end-to-end -p 16 -x $bowtie2Index -1 ${trim_IDs[0]} -2 ${trim_IDs[1]} | samtools view -bS - > ${sample_align_dir}/${sample_ID}.bam

echo "-------finish alignment for $sample_ID at $(date)--------" ;
}

total_sample=${#sample_IDs[@]}
n=0

for sample_ID in ${sample_IDs[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample samples---------------------- "

    task '$sample_ID' &

done

wait
echo "all done"


     
