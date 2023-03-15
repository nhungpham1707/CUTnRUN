#!/bin/bash
#SBATCH --job-name=bowtie2_all_samples
#SBATCH --output=all_samples_bowtie2.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=2
#SBATCH --mem=90G
#SBATCH --cpus-per-task=10
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl


# this script is to conduct trimming for all samples 

# Nhung, 09 03 2023

# Tool: bowtie2, ver 0.6.10

# Input: R1 and R2 outputs from trimming for each sample

# Output:
# - 1 bam file per sample stored in its respective folder

#########
# data dir



# sample_Ids was generated as text file from ls filename.txt from the data_dir
# remove       "SCC-bulkChIC-PMC-DRO-005"\ because it is running on test_bowtie2.sh to not overwrite
sample_IDs=( "1-bulkChIC-PMC-DRO-011" \
      "2-bulkChIC-PMC-DRO-013"\
      "3-bulkChIC-PMC-DRO-016"\
      "4-SCC-ChIC-PMC-DRO-L5"\
      "5-bulkChIC-PMC-DRO-012"\
      "6-SCC-bulkChIC-PMC-DRO-008"\
      "7-SCC-ChIC-PMC-DRO-LH"\
      "8-bulkChIC-PMC-DRO-013"\
      "9-SCC-ChIC-PMC-DRO-F1"\
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
   sleep $[($RANDOM % 10) + 1]s; echo "$1";
}

for str in ${sample_IDs[@]}; do
  #sleep $[ ( $RANDOM % 10 )  + 1 ]s; echo "$str"; &
  #n=$((n+1))
  #echo "-----------running $n out of $total_sample samples---------------------- "
  #echo $str &
  task "$str" &
done

wait"
#echo("all done")
