#!/usr/bin/env python3
import random

import pyrosetta as prs
from pyrosetta import *
import create_scfxn
from skopt.utils import use_named_args
from hyperparams import scfxn_ref15_space
from Bio.Align import substitution_matrices
from Bio import pairwise2

prs.init("-ex1 "
         "-ex2 "
         "-mute core.pack.pack_rotamers core.pack.task")  # no output from the design process

run = 0  # counter


@use_named_args(dimensions=scfxn_ref15_space)
def design_with_config(**config):
    global run
    run = run + 1
    print("RUN: ", run)

    print(config)
    ref15 = get_fa_scorefxn()  # REF15 score function with default weights for loss calulation
    scfxn = create_scfxn.creat_scfxn_from_config(config=config)  # optimization score Function
    pose = pose_from_pdb('benchmark/1K9P_A_relax_0001.pdb')  # easiest struct for optimizer picking
    native_pose = Pose()
    native_pose.assign(pose)
    resfile = "./design.resfile"
    with open(resfile, "w") as f:
        f.write("ALLAAxc \n")
        f.write("start\n")

    taskf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    taskf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    taskf.push_back(pyrosetta.rosetta.core.pack.task.operation.ReadResfile(resfile))
    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scfxn)
    packer.task_factory(taskf)
    taskf.create_task_and_apply_taskoperations(pose)
    packer.apply(pose)

    print("Optimized scfxn score: ", scfxn(pose))
    print('REF15 Score ', ref15(pose))
    # TODO: What defines our loss, for now use REF15
    bloss62 = substitution_matrices.load("BLOSUM62")
    similar = pairwise2.align.globaldx(pose.sequence(), native_pose.sequence(), bloss62, score_only=True)
    print('similarity ::', similar)
    return similar

