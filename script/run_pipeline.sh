#!/bin/bash
#SBATCH --job-name=module2
#SBATCH --output=module2.out
#SBATCH --time=96:0:0
#SBATCH --ntasks=1
#SBATCH --mem=90G
#SBATCH --cpus-per-task=8
#SBATCH --gres=tmpspace:30G
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=t.t.n.pham-3a@prinsesmaximacentrum.nl

# before running this script, change the Config.sh to fit your experiment. In addition, update variable in each module if needed.

set -a; source Config.sh  # set -a to export all variables to make it available for python and R (?)

# module 1. data processing: qc, trim, alignment, filtering, replicate correlation, bam2bw, merge replicates, peak calling for raw data, calculate frip 
. ./module-1.sh # replace condition group in peak calling with control and merge replicate steps

# module 2: data normalization: prepare data for normalization and conduct s3norm normalization

. ./module-2.sh

# module 3,5,6: are in a separete R script

# module 4: change variables in module-4.sh before running 
. ./module-4.sh
# module 7:
. ./module-7.sh