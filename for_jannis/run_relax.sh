#!/bin/bash

pdbtag=$1
num=$2
weights=$3

/dors/meilerlab/home/muellebk/rosetta/Rosetta_orbs_master/Rosetta/main/source/bin/relax.default.linuxgccrelease \
	@flags \
	@flags_$num \
	-s  /dors/meilerlab/home/muellebk/Projects/orbitals/wannier/dualoptE/decoys/xtal_refine/$pdbtag'_clean_0001.pdb' \
	-score:weights $weights \
	-set_weights cart_bonded 0.5 pro_close 0 yhh_planarity 0 \
	-nstruct 1 \
	-no_optH false -ignore_unrecognized_res -overwrite \
	-default_max_cycles 200 \
	-out:prefix opt_$num/ \
	-relax:script cartminpack.script \
	-relax:min_type lbfgs_armijo_nonmonotone \
	-extra_improper_file extra_improper_params.txt \
	-mute all
	#-crystal_refine \
	#-symmetry_definition /dors/meilerlab/home/muellebk/Projects/orbitals/wannier/dualoptE/decoys/xtal_refine/$pdbtag'.symm' \
