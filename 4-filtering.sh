#!/bin/bash
#SBATCH --job-name=filtering
#SBATCH --output=filtering.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=90G
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl


# this script is to sort and remove duplicate
# Nhung, 21 03 2023

# Tool: 
# - picard ver
# - samtools, ver 1.16.1, input flags:
# 	(https://wikis.utexas.edu/display/CoreNGSTools/Filtering+with+)
# 	(http://www.htslib.org/doc/samtools-view.html)
# 		-b to get output in bam file
# 		-q 20 skip alignment smaller than 20
# 		-F do not output alignment with any bit set in flag (1024 mean from pcr or optical dups)
# 		-h keep header
# 		-f only output alignment with all bits set in flag (2 = read map in proper pair https://www.samformat.info/sam-format-flag)

# Input data: 1 bam file for each sample (output from alignment w bowtie2)

# Output:
  # - 1 sorted bam file
  # - 1 mark and remove duplicate bam file
  # - 1 filter bam file (input for next step peakcalling)
  # - 1 index bam.bai file
  # - 1 metric text file

######### Define global variables ################

data_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/SCC_ChIC-PMC-DRO_plates_20210520_run1

# result dir
res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test

qc_dir=${res_dir}/qualityCheck
mkdir -p $qc_dir

trim_dir=${res_dir}/trimmed
mkdir -p $trim_dir

align_dir=${res_dir}/alignment
mkdir -p $align_dir

rm_dup_dir=$res_dir/rm_dup
mkdir -p $rm_dup_dir

peak_dir=${res_dir}/peakCalling/${sample_ID}
mkdir -p ${peak_dir}

peak_no_control_dir=${res_dir}/peakCalling_nocontrol
mkdir -p ${peak_no_control_dir}

# tool dir


bowtie2Index=/hpc/pmc_drost/SOURCES/Genomes/human/bowtie2/human_gencode37_hg38

picardTool=/hpc/pmc_drost/PROJECTS/swang/software/picard.jar
new_tmp_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/tmp # to solve the out of space with the temporary output from picard



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

########### Make task function ##############

task () {
  echo $sample_ID
    
  bam_ID=( $(find $align_dir/${sample_ID} -name "*.bam") )
  
  len=${#bam_ID[@]}
  echo $len
  sample_rmdup=$rm_dup_dir/${sample_ID}
  mkdir -p $sample_rmdup

  echo "start sorting and remove dup for $sample_ID at $(date)"
# sort
  java -Xmx2g -Djava.io.tmpdir=$new_tmp_dir -jar $picardTool SortSam \
	 -I ${bam_ID} \
	 -O ${sample_rmdup}/${sample_ID}_sorted.bam \
	 -VALIDATION_STRINGENCY LENIENT \
	 -SORT_ORDER coordinate \
	 -TMP_DIR $new_tmp_dir
	 echo "finsh sorting $sample_ID at $(date)"

# Mark and remove duplicates
  echo "start mark dup and remove for $sample_ID at $(date)"

  java -Xmx2g -Djava.io.tmpdir=$new_tmp_dir -jar $picardTool MarkDuplicates \
	 -I ${sample_rmdup}/${sample_ID}_sorted.bam  \
	 -O ${sample_rmdup}/${sample_ID}_mkrmdup.bam \
	 -VALIDATION_STRINGENCY LENIENT \
   	 -REMOVE_DUPLICATES true \
	 -M ${sample_rmdup}/${sample_ID}_marked_dup_metrics.txt \
	 -TMP_DIR $new_tmp_dir

  echo "finish mark and remove dup for $sample_ID at $(date)"

# Filter with samtools # 1024 for pcr or optical dups 

  echo "start sorting with samtools for $sample_ID at $(date)"

  samtools view -F 1024 -f 2 -q 20 -h -b ${sample_rmdup}/${sample_ID}_mkrmdup.bam -o ${sample_rmdup}/${sample_ID}_rmdup_filt.bam
# Index with samtools
# -@ n thread
  echo "start indexing for $sample_ID at $(date)"

  samtools index -@ 16 ${sample_rmdup}/${sample_ID}_rmdup_filt.bam ${sample_rmdup}/${sample_ID}_rmdup_filt.bam.bai

  echo "finish indexing for $sample_ID at $(date)" ;
}

############## run loop #################
n=0

for sample_ID in ${sample_IDs[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample samples---------------------- "

    task '$sample_ID' &

done

wait
echo "all done"














