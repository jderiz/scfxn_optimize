#!/bin/bash
# shellcheck disable=SC2206
#SBATCH --partition=galaxy-job
#SBATCH --job-name=CyTst
#SBATCH --output=out.log
#SBATCH --mem-per-cpu=2G
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-user=jannis.deriz@gmail.com
#SBATCH --time=48:00:00

## This script works for any number of nodes, Ray will find and manage all resources, if --exclusive
## Its Also possibe to only supply the number of cores to be used. IF this is the case one has to specifically pass the --cores(-c) flag to the script, else the script tries to distribute over all cores of each node that is being used.

#SBATCH --nodes=16
##SBATCH --ntasks=256
#SBATCH --ntasks-per-node=1 ## Ray can manage the Rsources
#SBATCH --cpus-per-task=24
#SBATCH --exclusive


##deprecated
###source /nfs/cluster/easybuild/software/Anaconda3/2020.02/etc/profile.d/conda.sh
# Load conda env
module load Anaconda3
source /software/all/Anaconda3/2020.02/etc/profile.d/conda.sh
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
  ray start --head --include-dashboard=true --node-ip-address="$ip" --port=$port --redis-password="$redis_password" --block &
sleep 30

worker_num=$((SLURM_JOB_NUM_NODES - 1)) #number of nodes other than the head node
for ((i = 1; i <= worker_num; i++)); do
  node_i=${nodes_array[$i]}
  echo "STARTING WORKER $i at $node_i"
  srun --nodes=1 --ntasks=1 -w "$node_i" ray start --address "$ip_head" --redis-password="$redis_password" --block &
  sleep 3
done
echo "FORWARD DASHBOARD PORT 8265 TO LOCAL MACHIN FOR DASHBOARD"

# ===== Call your code below =====
prot=$1
target=$2
loss=$3
python bench.py -e RF  -l "$loss" -rpc 6 -cycles 3 -evals 20 -pdb "$prot" -target "$target" -id optimize_"$prot"_"$loss" -cooldown -r_pw "$redis_password"
