#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=2G
#SBATCH --partition=clara-job
#SBATCH --mail-type=end
#SBATCH --mail-user=jannis.deriz@gmail.com
#SBATCH -n 1
#SBATCH -c 64
#SBATCH --job-name=simi_gp


module load Anaconda3
source /nfs/cluster/easybuild/software/Anaconda3/2020.02/etc/profile.d/conda.sh

conda activate scfxn

python -u bench.py GP > py_out/GP.out

