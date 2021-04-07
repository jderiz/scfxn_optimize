#!/bin/sh

pdbtag=$1
num=$2
weights=$3

rm -f $outdir/$pdbtag*.sc.min

mkdir -p opt_$num/$pdbtag'_min'

/dors/meilerlab/home/muellebk/rosetta/Rosetta_orbs_master/Rosetta/main/source/bin/min_test.default.linuxgccrelease \
	@flags \
	@flags_$num \
	-in:file:s  /dors/meilerlab/home/muellebk/Projects/orbitals/wannier/dualoptE/decoys/monomer_pdbs_out_trim/$pdbtag*.pdb.gz \
	-score:weights $weights \
	-set_weights cart_bonded 0.5 pro_close 0 yhh_planarity 0 \
	-min:minimizer lbfgs_armijo_nonmonotone \
	-min::pack \
	-default_max_cycles 50 \
	-min:cartesian \
	-out:path:score opt_$num/ \
	-out:file:scorefile opt_$num/$pdbtag.sc.min.no_rmsd \
	-out:path:pdb opt_$num/$pdbtag'_min'/ \
	-mute all

/dors/meilerlab/home/muellebk/rosetta/Rosetta_orbs_master/Rosetta/main/source/bin/rosetta_scripts.default.linuxgccrelease \
	@flags \
	@flags_$num \
	-in:file:s opt_$num/$pdbtag'_min'/*.pdb \
	-parser:protocol min_test.xml \
	-parser:script_vars wts=$weights \
	-score:weights $weights \
	-set_weights cart_bonded 0.5 pro_close 0 yhh_planarity 0 \
	-out:file:score_only opt_$num/$pdbtag.sc.min \
	-in:file:native  /dors/meilerlab/home/muellebk/Projects/orbitals/wannier/dualoptE/decoys/monomer/$pdbtag'_clean.pdb' \
	-mute all
