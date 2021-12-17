#!/bin/bash
for prot in "6q21A" "1avsA" "1lfaA"
do
sbatch runc.sh $prot rmsd
sbatch runc.sh $prot score
done
