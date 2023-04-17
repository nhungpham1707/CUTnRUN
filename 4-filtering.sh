#!/bin/bash

# this script is to sort and remove duplicate
# Nhung, 21 03 2023

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

# task () {
#   echo $sample_ID
    
#   bam_ID=( $(find $align_dir/${sample_ID} -name "*.bam") )
  
#   len=${#bam_ID[@]}
#   echo $len
#   sample_rmdup=$rm_dup_dir/${sample_ID}
#   mkdir -p $sample_rmdup

#   echo "start sorting and remove dup for $sample_ID at $(date)"
# # sort
#   java -jar $picardTool SortSam \
# 	 -I ${bam_ID} \
# 	 -O ${sample_rmdup}/${sample_ID}_sorted.bam \
# 	 -VALIDATION_STRINGENCY LENIENT \
# 	 -SORT_ORDER coordinate 
	 
#   echo "finsh sorting $sample_ID at $(date)"

# # Mark and remove duplicates
#   echo "start mark dup and remove for $sample_ID at $(date)"

#   java -jar $picardTool MarkDuplicates \
# 	 -I ${sample_rmdup}/${sample_ID}_sorted.bam  \
# 	 -O ${sample_rmdup}/${sample_ID}_mkrmdup.bam \
# 	 -VALIDATION_STRINGENCY LENIENT \
#    	 -REMOVE_DUPLICATES true \
# 	 -M ${sample_rmdup}/${sample_ID}_marked_dup_metrics.txt 

#   echo "finish mark and remove dup for $sample_ID at $(date)"

# # Filter with samtools # 1024 for pcr or optical dups 

#   echo "start sorting with samtools for $sample_ID at $(date)"

#   samtools view -F 1024 -f 2 -q 20 -h -b ${sample_rmdup}/${sample_ID}_mkrmdup.bam -o ${sample_rmdup}/${sample_ID}_rmdup_filt.bam
# # Index with samtools
# # -@ n thread
#   echo "start indexing for $sample_ID at $(date)"

#   samtools index -@ 16 ${sample_rmdup}/${sample_ID}_rmdup_filt.bam ${sample_rmdup}/${sample_ID}_rmdup_filt.bam.bai

#   echo "finish indexing for $sample_ID at $(date)" ;
# }
############## run loop #################
n=0

for sample_ID in ${sample_IDs[@]}; do
    n=$((n+1))
    echo "-----------running $n out of $total_sample samples---------------------- "

    task '$sample_ID' &

done

wait

# generate report 
multiqc ${rm_dup_dir} -n rm_dup_report -o ${rm_dup_dir}
echo "all done filtering"














