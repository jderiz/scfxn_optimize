#!/bin/bash
declare -A bench_pairs
bench_pairs=( [../benchmark/allosteric/pre_relaxed/no_ligand_relaxed_starting_1k9kA.pdb]=../benchmark/allosteric/pre_relaxed/relaxed_target_1K9P.pdb)
#              [../benchmark/allosteric/pre_relaxed/no_ligand_relaxed_starting_3s0bA.pdb]=../benchmark/allosteric/pre_relaxed/relaxed_target_3S0A.pdb
#              [../benchmark/allosteric/pre_relaxed/no_ligand_relaxed_starting_1f4vA.pdb]=../benchmark/allosteric/pre_relaxed/relaxed_target_3CHY.pdb
#              [../benchmark/allosteric/pre_relaxed/no_ligand_relaxed_starting_3zjaA.pdb]=../benchmark/allosteric/pre_relaxed/relaxed_target_3ZK0.pdb
#              [../benchmark/allosteric/pre_relaxed/no_ligand_relaxed_starting_6q21A.pdb]=../benchmark/allosteric/pre_relaxed/relaxed_target_4Q21.pdb
#              [../benchmark/allosteric/pre_relaxed/no_ligand_relaxed_starting_1lfaA.pdb]=../benchmark/allosteric/pre_relaxed/relaxed_target_1MQ9.pdb
#              [../benchmark/allosteric/pre_relaxed/no_ligand_relaxed_starting_1d5wA.pdb]=../benchmark/allosteric/pre_relaxed/relaxed_target_1D5B.pdb
# iterate over map
for prot in ${!bench_pairs[@]}
do
echo "run $prot to ${bench_pairs[$prot]}"
sbatch runc.sh ${prot} ${bench_pairs[$prot]} score
done
