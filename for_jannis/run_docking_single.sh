#!/bin/sh

pdbtag=$1
num=$2
weights=$3

rm -f $outdir/$pdbtag*.ddg.out

/dors/meilerlab/home/muellebk/rosetta/Rosetta_orbs_master/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
	-parser:protocol ddg.xml \
	-parser:script_vars wts=$weights \
	@flags \
	@flags_$num \
	-in:file:s /dors/meilerlab/home/muellebk/Projects/orbitals/wannier/dualoptE/decoys/docking_pdbs/$pdbtag*.pdb.gz \
	-score:weights $weights \
	-in:file:native  /dors/meilerlab/home/muellebk/Projects/orbitals/wannier/dualoptE/decoys/docking/$pdbtag'_bound_native.pdb' \
	-out:file:score_only opt_$num/$pdbtag.ddg.out  \
	-silent_read_through_errors \
	-mute all
