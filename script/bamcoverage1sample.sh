#!/bin/bash
#SBATCH --job-name=bamcoverage
#SBATCH --output=bamcoverage_2sample.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --cpus-per-task=16
#SBATCH --gres=tmpspace:30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl

# nhung 25 04 2023
# test bamcoverage bedgraph file for one sample. if the bedgraph is readable (by igv or head ss.bedgraph) with the outputformat define
res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test
rm_dup_dir=$res_dir/rm_dup
bincount_dir=${res_dir}/bincount_bamcoverage_2sample
mkdir -p $bincount_dir
bamCoverage_dir=/hpc/pmc_drost/nhung/anaconda3/envs/cutnrun_trimgalore/bin/bamCoverage # to merge 
effectiveGenomeSize=2913022398
bin_size=200
# sample_ID=bulkChIC-PMC-DRO-011
sample_IDs=( "bulkChIC-PMC-DRO-011" "SCC-ChIC-PMC-DRO-T1")
for sample_ID in ${sample_IDs[@]}
do
$bamCoverage_dir -b ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam \
        --binSize $bin_size \
        --effectiveGenomeSize $effectiveGenomeSize \
        --numberOfProcessors 16 \
        --outFileFormat bedgraph \
        -o ${bincount_dir}/${sample_ID}_binsize_${bin_size}.bedgraph
done