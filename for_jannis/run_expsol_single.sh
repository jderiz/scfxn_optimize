#!/bin/bash

pdbtag=$1
num=$2
weights=$3

~/rosetta/Rosetta_vanilla/Rosetta/main/source/bin/fasol_perres.default.linuxgccrelease \
	@flags \
	@flags_$num \
	-score:weights $weights \
	-s /dors/meilerlab/home/muellebk/Projects/orbitals/wannier/dualoptE/decoys/fasol/$pdbtag.pdb \
	-mute all -unmute fasol_refit \
	-overwrite -nooutput \
	-ignore_unrecognized_res \
	-crystal_refine \
	-flip_HNQ -no_optH false > opt_$num/$pdbtag.EXPSOL
