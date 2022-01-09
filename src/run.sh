#!/bin/bash
for prot in "1d5wA" "6q21A"  "1lfaA" "1f4vA" "1k9kA" "3zjaA" "3s0bA"
do
#sbatch runc.sh $prot score
sbatch run_FRB.sh $prot
done
