#!/bin/bash
#SBATCH --time=32:00:00
#SBATCH --mem=10G
#SBATCH --mail-type=FAIL
#SBATCH --partition=clara-job
#SBATCH --mail-user=jannis.deriz@gmail.com
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --job-name=simi_forest

module load Anaconda3
source /nfs/cluster/easybuild/software/Anaconda3/2020.02/etc/profile.d/conda.sh

conda activate scfxn

python -u bench.py RF > py_out/RF.out

