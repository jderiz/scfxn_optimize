#!/bin/bash

pdbtag=$1
partnum=$2
num=$3
weights=$4

/dors/meilerlab/home/muellebk/rosetta/Rosetta_orbs_master/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
	@flags \
	@flags_$num \
	-s /dors/meilerlab/home/muellebk/Projects/orbitals/wannier/dualoptE/decoys/rotrecov/core/$pdbtag'_clean.pdb' \
	-score:weights $weights \
	-out::nooutput \
	-mapfile /dors/meilerlab/home/muellebk/Projects/orbitals/wannier/dualoptE/decoys/rotrecov/core/$pdbtag.ccp4 \
	-parser:protocol features.xml \
	-parser:script_vars outdir=opt_$num pdb=$pdbtag resfile=/dors/meilerlab/home/muellebk/Projects/orbitals/wannier/dualoptE/decoys/rotrecov/core/$pdbtag.exclude_multiple_exclude_xtal_exclude_density_exclude_bfactors.resfile database_partition=$partnum type=core \
	-mute all #-ignore_unrecognized_res

