#!/bin/sh

pdbtag=$1
num=$2
weights=$3

rm -f $outdir/$pdbtag*.sc
/dors/meilerlab/home/muellebk/rosetta/Rosetta_orbs_master/Rosetta/main/source/bin/score_jd2.default.linuxgccrelease \
	@flags \
	@flags_$num \
	-in:file:s /dors/meilerlab/home/muellebk/Projects/orbitals/wannier/dualoptE/decoys/monomer_pdbs_out/$pdbtag*.pdb.gz \
	-out:file:scorefile opt_$num/$pdbtag.sc.static \
	-score:weights $weights \
	-set_weights cart_bonded 0.5 pro_close 0 yhh_planarity 0 \
	-in:file:native /dors/meilerlab/home/muellebk/Projects/orbitals/wannier/dualoptE/decoys/monomer/$pdbtag'_clean.pdb' \
	-mute all
