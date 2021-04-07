#!/bin/sh

# core RR
bash scripts/merge.sh $1/features.core.db3 $1/*core.db3_*
LD_PRELOAD=/lib64/libsqlite3.so.0 /usr/bin/sqlite3 $1/features.core.db3  < scripts/average_rotamer_recovery_by_amino_acid.sql > $1/rr_core_by_amino_acid.csv
LD_PRELOAD=/lib64/libsqlite3.so.0 /usr/bin/sqlite3 $1/features.core.db3 < scripts/average_rotamer_recovery.sql > $1/rr_core_result

# int RR
bash scripts/merge.sh $1/features.int.db3 $1/*int.db3_*
LD_PRELOAD=/lib64/libsqlite3.so.0 /usr/bin/sqlite3 $1/features.int.db3  < scripts/average_rotamer_recovery_by_amino_acid.sql > $1/rr_interface_by_amino_acid.csv
LD_PRELOAD=/lib64/libsqlite3.so.0 /usr/bin/sqlite3 $1/features.int.db3 < scripts/average_rotamer_recovery.sql > $1/rr_interface_result

# rescore / score_decoy
ls $1/*.sc.static > temp.scores
./calc1dboltzmann.pl score temp.scores > $1/score_decoy_result
rm $1/*.sc.static

# mintest / score_min
ls $1/*.sc.min > temp.scoresmin
./calc1dboltzmann.pl total_score temp.scoresmin > $1/score_min_result
rm -r $1/*_min/
rm $1/*.sc.min.no_rmsd
rm $1/*.sc.min

# ddgs / docking_single
ls $1/*.ddg.out > temp.ddg
./calc1dboltzmann.pl ddg temp.ddg > $1/score_ddg_result
rm $1/*.ddg.out

# xtal grads
#cat $1/*.grads | awk '{ sum += $1; n++ } END { if (n > 0) print sum/n, n; }' > $1/xtalgrad_result
#rm -rf _ros* $1/*.grads

# seqrecov / ssm_single
ln -s /home/muellebk/Projects/orbitals/wannier/dualoptE/decoys/seqrecov/*COUNTS $1/
./fitref `ls $1/*.ENERGIES | sed 's/\.ENERGIES//'` > $1/seqrecov_result
rm $1/*COUNTS $1/*ENERGIES

# expsol
#cat $1/*.EXPSOL > $1/ALL.EXPSOL
#./exp_desolvation.pl $1/ALL.EXPSOL > $1/expsol_result
#rm $1/*EXPSOL

# distrs / relax
./distdstr.py $1 > $1/distr_result
rm $1/*_0001.pdb

# ligand docking
sh ligdockanalysis.sh $1/ligdock/
tail -n+2 $1/ligdock/all_pnear.txt | awk '{ sum += $2 } END { if (NR > 0) print sum / NR }' > $1/ligdock_result
#rm -r $1/ligdock/

# final cleanup
rm temp.ddg temp.scoresmin temp.scores
rm rosetta_array_job_slurm_*.out
