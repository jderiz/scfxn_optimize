#!/bin/bash

pdbtag=$1
partnum=$2
num=$3
weights=$4

/dors/meilerlab/home/muellebk/rosetta/Rosetta_orbs_master/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
	@flags \
	@flags_$num \
	-s /dors/meilerlab/home/muellebk/Projects/orbitals/wannier/dualoptE/decoys/rotrecov/interface/$pdbtag'_clean.pdb' \
	-score:weights $weights \
	-mapfile /dors/meilerlab/home/muellebk/Projects/orbitals/wannier/dualoptE/decoys/rotrecov/interface/$pdbtag.map \
	-out::nooutput \
	-parser:protocol features.xml \
	-parser:script_vars outdir=opt_$num pdb=$pdbtag resfile=/dors/meilerlab/home/muellebk/Projects/orbitals/wannier/dualoptE/decoys/rotrecov/interface/$pdbtag"_cleanwat_0001_eval.resfile" database_partition=$partnum type=int \
	-mute all #-ignore_unrecognized_res

