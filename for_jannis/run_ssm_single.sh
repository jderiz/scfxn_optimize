#!/bin/bash

pdbtag=$1
num=$2
weights=$3

/dors/meilerlab/home/muellebk/rosetta/Rosetta_orbs_master/Rosetta/main/source/bin/calc_ssm_energies.default.linuxgccrelease \
	@flags \
	@flags_$num \
	-s /dors/meilerlab/home/muellebk/Projects/orbitals/wannier/dualoptE/decoys/seqrecov/$pdbtag'_0001.pdb' \
	-score:weights $weights \
	-out::nooutput \
	-mute all -unmute calc_ssm_energies > opt_$num/$pdbtag.ENERGIES
