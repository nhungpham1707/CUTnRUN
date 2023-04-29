#!/bin/bash
#SBATCH --job-name=bamcoverage
#SBATCH --output=bamcoverage_rerun_1stsample.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --cpus-per-task=16
#SBATCH --gres=tmpspace:30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl

# nhung 25 04 2023
# test bamcoverage bedgraph file for one sample. if the bedgraph is readable (by igv or head ss.bedgraph) with the outputformat define
res_dir=/hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test
rm_dup_dir=$res_dir/rm_dup
bincount_dir=${res_dir}/bincount_bamcoverage_fixbinsize
mkdir -p $bincount_dir
bamCoverage_dir=/hpc/pmc_drost/nhung/anaconda3/envs/cutnrun_trimgalore/bin/bamCoverage # to merge 
effectiveGenomeSize=2913022398
bin_size=200
sample_IDs=bulkChIC-PMC-DRO-011
# sample_IDs=( "bulkChIC-PMC-DRO-011" "SCC-ChIC-PMC-DRO-T1")
# sample_IDs=( "bulkChIC-PMC-DRO-011" \
#             "bulkChIC-PMC-DRO-012"\
#             "bulkChIC-PMC-DRO-013"\
#             "bulkChIC-PMC-DRO-014"\
#             "bulkChIC-PMC-DRO-015"\
#             "bulkChIC-PMC-DRO-016"\
#             "SCC-bulkChIC-PMC-DRO-002"\
#             "SCC-bulkChIC-PMC-DRO-005"\
#             "SCC-bulkChIC-PMC-DRO-008"\
#             "SCC-ChIC-PMC-DRO-L5"\
#             "SCC-ChIC-PMC-DRO-LH"\
#             "SCC-ChIC-PMC-DRO-F1"\
#             "SCC-ChIC-PMC-DRO-F5"\
#             "SCC-ChIC-PMC-DRO-FH"\
#             "SCC-ChIC-PMC-DRO-T1"\
#             "SCC-ChIC-PMC-DRO-T5"\
#             "SCC-ChIC-PMC-DRO-TH"\ 
#             "SCC-ChIC-PMC-DRO-L1"\ 
#             "SCC-bulkChIC-PMC-DRO-020"\
#             "SCC-bulkChIC-PMC-DRO-021"\
#             "SCC-bulkChIC-PMC-DRO-022"\
#             "SCC-bulkChIC-PMC-DRO-023"\
#             "SCC-bulkChIC-PMC-DRO-024"\
#             "SCC-bulkChIC-PMC-DRO-025")

for sample_ID in ${sample_IDs[@]}
do
$bamCoverage_dir -b ${rm_dup_dir}/${sample_ID}/${sample_ID}_rmdup_filt.bam \
        --binSize $bin_size \
        --effectiveGenomeSize $effectiveGenomeSize \
        --numberOfProcessors 16 \
        --outFileFormat bedgraph \
        -o ${bincount_dir}/${sample_ID}_binsize_${bin_size}.bedgraph
done