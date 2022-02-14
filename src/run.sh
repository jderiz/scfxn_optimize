#!/bin/bash
declare -A bench_pairs
bench_pairs=([1k9k_A.pdb]=1k9p_A.pdb
              [3s0b_A.pdb]=3s0a_A.pdb
              [1f4v_A.pdb]=3chy_A.pdb
              [3zja_A.pdb]=3zk0_A.pdb
              [1d5w_A.pdb]=1dbw_A.pdb)

# iterate over map of bound to unbound
for prot_path in ${!bench_pairs[@]}
do
  # get prot name
prot=${prot_path: -10:6}
echo "PDB:" ${prot}
echo "run $prot_path to ${bench_pairs[$prot_path]}"
# run sbatch with correct paths
# get print availiable cpus
cluster="clara-job"
avail_cpus=$(sinfo -p "${cluster}" -o "%C" | tail -n1 | grep -Eo "[0-9]+" | head -n2 | tail -n1)
echo "IDLE CPUS ON CLUSTER ${cluster}: ${avail_cpus}"
sbatch runc.sh ${prot} ../benchmark/allosteric/pre_processed/${prot_path} ../benchmark/allosteric/pre_processed/${bench_pairs[$prot_path]} ca_rmsd
done
