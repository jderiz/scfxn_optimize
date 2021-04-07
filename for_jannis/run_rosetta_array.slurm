#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH --ntasks=1
#SBATCH --time=2:00:00
#SBATCH --mem=2G
#SBATCH --array=1-288   # go from 1 to the number of files in JOB_LIST
#SBATCH --output=rosetta_array_job_slurm_%A_%a.out

export LD_LIBRARY_PATH=/dors/meilerlab/apps/Linux2/x86_64/lib64:/dors/meilerlab/apps/Linux2/x86_64/lib:$LD_LIBRARY_PATH

arrayfile=`cat ${JOB_LIST} | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR == line) print $0}'`

sh $arrayfile ${NUM} ${WEIGHTS}
