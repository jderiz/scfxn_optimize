#!/bin/sh

pdbtag=$1
num=$2
weights=$3

phenix.rosetta.run_phenix_interface \
	~/Rosetta/main/source/bin/rosetta_scripts.python.linuxgccrelease \
	-parser:protocol refine.xml \
	-parser:script_vars symmdef=../decoys/xtal_refine/$pdbtag.symm wts=$weights outfile=opt_$num/$pdbtag.grads  \
	@flags \
	@flags_$num \
	-s ../decoys/xtal_refine/$pdbtag'_clean_0001.pdb' \
	-cryst::mtzfile ../decoys/xtal_refine/$pdbtag'-sf.mtz' \
	-crystal_refine \
	-nstruct 1 \
	-no_optH false -ignore_unrecognized_res -overwrite \
	-renumber_pdb false \
	-out::nooutput \
	-mute all #-unmute protocols.cryst.cryst_movers &> opt_$num/$pdbtag.grad_detail
