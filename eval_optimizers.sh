#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --mem=1G
#SBATCH --mail-type=end
#SBATCH --mail-user=jannis.deriz@gmail.com
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --job-name=eval_optimizers


module load Anaconda3
source /nfs/cluster/easybuild/software/Anaconda3/2020.02/etc/profile.d/conda.sh

conda activate scfxn

python bench.py forest

