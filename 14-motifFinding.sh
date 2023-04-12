#!/bin/bash

# This script is to find motif from macs2 narrow peak 
# more reference for flags or output interpretation here http://homer.ucsd.edu/homer/ngs/peakMotifs.html
# Nhung 22 03 2023


task () {

out_dir=${motif_dir}/${sample_ID}
mkdir -p ${out_dir}

sample_dir=( $(find ${peak_dir}/${sample_ID}/narrow -name "*paired_control_summits.bed") ) 
echo "$sample_dir"
echo "----------------start-motif-finding-for $sample_ID at $(date)----------"

#$findMotif_dir $sample_dir hg38 $_out_dir
# $findMotif_dir $sample_dir hg38 peakAnalysis -size 200 -len 8 $_out_dir

# findMotifsGenome.pl $sample_dir $fasta_genome_dir peakAnalysis -size 200 -len 8 $_out_dir

findMotifsGenome.pl $sample_dir $fasta_genome_dir $out_dir -size 200 -len 8 

#findMotifsGenome.pl ERalpha.peaks hg18 MotifOutputDirectory/ -find motif1.motif > outputfile.txt
echo "----------------finish-motif-finding-for $sample_ID at $(date)----------" ;
}
# cannot read file it seems. check solution here https://www.biostars.org/p/269709/

total_sample=${#motif_samples[@]}
n=0
for sample_ID in $motif_samples; do

n=$((n+1))
  echo "-----------running $n out of $total_sample  samples---------------"

  task "$sample_ID" &

done

wait

echo "all done"
