#!/bin/bash
#SBATCH --job-name=func_paral_bowtiee
#SBATCH --output=func_paral_all_sample_bowtie2_bam_file.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=3
#SBATCH --mem=40G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl


# this script is to conduct alignment for all samples 

# Nhung, 15 03 2023

# Tool: bowtie2, ver 0.6.10

# Input: R1 and R2 outputs from trimming for each sample

# Output:
# - 1 bam file per sample stored in its respective folder

#########
# data dir
data_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/SCC_ChIC-PMC-DRO_plates_20210520_run1
res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test
mkdir -p ${res_dir}/alignment
align_dir=${res_dir}/alignment

genomelib=/hpc/pmc_gen/references/RNA-Seq/starfusion/GRCh38_gencode_v31_CTAT_lib_Oct012019.plug-n-play/ctat_genome_lib_build_dir/
gtf=/hpc/pmc_drost/SOURCES/Genomes/human/gencode37_GRCh38_annotation.gtf
bowtie2Index=/hpc/pmc_drost/SOURCES/Genomes/human/bowtie2/human_gencode37_hg38


# sample_Ids was generated as text file from ls filename.txt from the data_dir
sample_IDs=( "bulkChIC-PMC-DRO-011"\
	     "bulkChIC-PMC-DRO-013"\
      "bulkChIC-PMC-DRO-016"\
      "SCC-bulkChIC-PMC-DRO-005"\
      "SCC-ChIC-PMC-DRO-L5"\
      "bulkChIC-PMC-DRO-012"\
      "SCC-bulkChIC-PMC-DRO-008"\
      "SCC-ChIC-PMC-DRO-LH"\
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
    
total_sample=${#sample_IDs[@]}
n=0

task(){
	trim_IDs=( $(find ${res_dir}/trimmed/${sample_ID} -name "*.fq.gz") )
  #declare -p fastq_IDs

  len=${#trim_IDs[@]}
  echo $len
  
  mkdir -p ${align_dir}/${sample_ID}
  sample_align_dir=${align_dir}/${sample_ID}
  
  echo "-----start alignment for $sample_ID at $(date)-------"

##Alignment to both target and spike-in genomes (Bowtie 2)
#cd $sample_align_dir # go to where the res will be saved

bowtie2 --end-to-end -p 16 -x $bowtie2Index -1 ${trim_IDs[0]} -2 ${trim_IDs[1]} | samtools view -bS - > ${sample_align_dir}/${sample_ID}.bam

  echo "-------finish alignment for $sample_ID at $(date)--------"; 	}


for str in ${sample_IDs[@]}; do
  n=$((n+1))
  echo "-----------running $n out of $total_sample samples---------------------- "
  sample_ID=$str
  task "$sample_ID" &	
 
done

wait

echo "all done"

