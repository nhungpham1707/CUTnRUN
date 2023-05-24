#!/bin/bash
#SBATCH --job-name=bw2bed
#SBATCH --output=convert_bigwig_to_bed.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --cpus-per-task=8
#SBATCH --gres=tmpspace:30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl

cd /hpc/pmc_drost/PROJECTS/swang/CUT_RUN/nhung_test/merged_bigwig

# bigWigToWig fusion_merged.bw fusion_merged.wig
wig2bed < fusion_merged.wig > fusion_merged.bed

# bigWigToWig tfe3_merged.bw tfe3_merged.wig
wig2bed < tfe3_merged.wig > tfe3_merged.bed

# bigWigToWig luciferase_merged.bw luciferase_merged.wig
wig2bed < luciferase_merged.wig > luciferase_merged.bed
