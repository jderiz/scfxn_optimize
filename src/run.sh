#!/bin/bash
declare -A bench_pairs
bench_pairs=( [1k9k_A.pdb]=1k9p_A.pdb)
#              [3s0b_A.pdb]=3s0a_A.pdb
#              [1f4v_A.pdb]=3chy_A.pdb
#              [3zja_A.pdb]=3zk0_A.pdb
#              [1d5w_A.pdb]=1dbw_A.pdb

# iterate over map of bound to unbound
for prot_path in ${!bench_pairs[@]}
do
  # get prot name
prot=${prot_path: -9:5}
echo "PDB:" ${prot}
echo "run $prot_path to ${bench_pairs[$prot_path]}"
# run sbatch with correct paths
sbatch runc.sh ${prot} ${prot_path} ${bench_pairs[$prot_path]} rmsd
done
