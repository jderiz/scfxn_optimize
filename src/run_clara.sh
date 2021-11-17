#!/bin/bash
# shellcheck disable=SC2206
#SBATCH --partition=clara-job
#SBATCH --job-name=AlostK9
#SBATCH --output=clara.log
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=jannis.deriz@gmail.com
#SBATCH --time=48:00:00
## This script works for any number of nodes, Ray will find and manage all resources
## Its Also possibe to only supply the number of cores to be used. IF this is the case one has to specifically pass the --cores(-c) flag to the script, else the script tries to distribute over all cores of each node that is being used.
#SBATCH --ntasks=256
## Give all resources to a single Ray task, ray can manage the resources internally
#SBATCH --gpus-per-task=0

# Load modules or your own conda environment here
# module load pytorch/v1.4.0-gpu
# conda activate ${CONDA_ENV}
module load Anaconda3
source /nfs/cluster/easybuild/software/Anaconda3/2020.02/etc/profile.d/conda.sh
conda activate scfxn

# ===== DO NOT CHANGE THINGS HERE UNLESS YOU KNOW WHAT YOU ARE DOING =====
# This script is a modification to the implementation suggest by gregSchwartz18 here:
# https://github.com/ray-project/ray/issues/826#issuecomment-522116599
redis_password=$(uuidgen)
export redis_password
printf $SLURM_JOB_NODELIST

nodes=$(scontrol show hostnames "$SLURM_JOB_NODELIST") # Getting the node names
nodes_array=($nodes)
node_1=${nodes_array[0]}
ip=$(srun --nodes=1 --ntasks=1 -w "$node_1" hostname --ip-address) # making redis-address

# if we detect a space character in the head node IP, we'll
# convert it to an ipv4 address. This step is optional.
if [[ "$ip" == *" "* ]]; then
  IFS=' ' read -ra ADDR <<< "$ip"
  if [[ ${#ADDR[0]} -gt 16 ]]; then
    ip=${ADDR[1]}
  else
    ip=${ADDR[0]}
  fi
  echo "IPV6 address detected. We split the IPV4 address as $ip"
fi

port=6379
ip_head=$ip:$port
export ip_head
echo "IP Head: $ip_head"

echo "STARTING HEAD at $node_1"
srun --nodes=1 --ntasks=1 -w "$node_1" \
  ray start --head --include-dashboard true --node-ip-address="$ip" --port=$port --redis-password="$redis_password" --block &
sleep 30

worker_num=$((SLURM_JOB_NUM_NODES - 1)) #number of nodes other than the head node
for ((i = 1; i <= worker_num; i++)); do
  node_i=${nodes_array[$i]}
  echo "STARTING WORKER $i at $node_i"
  srun --nodes=1 --ntasks=1 -w "$node_i" ray start --address "$ip_head" --redis-password="$redis_password" --block &
  sleep 5
done
echo "FORWARD DASHBOARD PORT 8265"
ssh -N -L 8265:localhost:8265 $node_1
# ===== Call your code below =====
python bench.py -l score -c 256 -rpc 6 -evals 600 -pdb 1f4vA -id allost_1f4vA -cooldown -r_pw "$redis_password"
