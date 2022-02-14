#!/bin/bash
# shellcheck disable=SC2206
#SBATCH --partition=clara-cpu
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
##NodeName=galaxy120 Arch=x86_64 CoresPerSocket=6
##   CPUAlloc=24 CPUTot=24 CPULoad=12.02
##      AvailableFeatures=(null)
##         ActiveFeatures=(null)
##            Gres=(null)
##               NodeAddr=galaxy120 NodeHostName=galaxy120 Version=21.08.5
##                  OS=Linux 3.10.0-1160.49.1.el7.x86_64 #1 SMP Tue Nov 30 15:51:32 UTC 2021
##                     RealMemory=125952 AllocMem=36864 FreeMem=84876 Sockets=2 Boards=1
##                        State=ALLOCATED ThreadsPerCore=2 TmpDisk=0 Weight=1 Owner=N/A MCS_label=N/A
##                           Partitions=galaxy-job
##                              BootTime=2022-01-12T07:23:44 SlurmdStartTime=2022-01-12T07:26:24
##                                 LastBusyTime=2022-01-25T03:37:42
##                                    CfgTRES=cpu=24,mem=123G,billing=24
##                                       AllocTRES=cpu=24,mem=36G
##                                          CapWatts=n/a
##                                             CurrentWatts=0 AveWatts=0
##                                                ExtSensorsJoules=n/s ExtSensorsWatts=0 ExtSensorsTemp=n/s
##                                                
##SBATCH --nodes=16
#SBATCH --ntasks=192
##SBATCH --ntasks-per-node=1 ## Ray can manage the Rsources
##SBATCH --cpus-per-task=24
##SBATCH --exclusive

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
printf "$SLURM_JOB_NODELIST"


nodes=$(scontrol show hostnames "$SLURM_JOB_NODELIST") # Getting the node names
nodes_array=($nodes) #make array from names
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

# make array from comma separated list
cpus_per_node=(${SLURM_JOB_CPUS_PER_NODE//,/ }) 
# parse (xN) type

echo "CPUS_PER_NODE: ${cpus_per_node}"
echo "STARTING HEAD at $node_1"
srun --nodes=1 --ntasks=1 -w "$node_1" \
  ray start --head --num-cpus "${cpus_per_node[0]}" --include-dashboard=true --node-ip-address="$ip" --port=$port --redis-password="$redis_password" --block &
sleep 3

worker_num=$((SLURM_JOB_NUM_NODES - 1)) #number of nodes other than the head node
for ((i = 1; i <= worker_num; i++)); do
  node_i=${nodes_array[$i]}
  echo "STARTING WORKER $i at $node_i"
  srun --nodes=1 --ntasks=1 -w "$node_i" ray start --num-cpus "${cpus_per_node[$i]}" --address "$ip_head" --redis-password="$redis_password" --block &
  sleep 3
done
echo "FORWARD DASHBOARD PORT 8265 TO LOCAL MACHIN FOR DASHBOARD"

# ===== Call your code below =====
prot=$1
prot_path=$2
target=$3
loss=$4
python bench.py -config ""  -l "$loss" -c 40 -evals 40 -pdb "$prot_path" -target "$target" -id FRB_40_${prot} -r_pw "$redis_password"
