#!/bin/bash
for prot in "6q21A"  "1lfaA" "1f4vA" "1k9kA" "1d5WA" "3zjaA" "3s0bA"
do
sbatch runc.sh $prot score
#sbatch runc.sh $prot rmsd
done
