#!/usr/bin/env python3
import random, os

import pyrosetta as prs
from Bio import pairwise2
from Bio.Align import substitution_matrices
from pyrosetta import (Pose, PyJobDistributor, get_fa_scorefxn, init,
                       pose_from_pdb)
from skopt.utils import use_named_args

import create_scfxn
from hyperparams import scfxn_ref15_space

prs.init(
    "-ex1 " "-ex2 " "-mute core.pack.pack_rotamers core.pack.task"
)  # no output from the design process

run = 0  # counter


@use_named_args(dimensions=scfxn_ref15_space)
def design_with_config(**config):
    global run
    run = run + 1
    print("RUN: ", run)

    print(config)
    ref15 = (
        get_fa_scorefxn()
    )  # REF15 score function with default weights for loss calulation
    scfxn = create_scfxn.creat_scfxn_from_config(
        config=config
    )  # optimization score Function
    pose = pose_from_pdb(
        "benchmark/1K9P_A_relax_0001.pdb"
    )  # easiest struct for estiamtor picking
    # TODO: Add random struct each run once decided on estimator
    # pdbs = []

    # for pdb in os.listdir("benchmark"):
    #     if pdb.endswith(".pdb"):
    #         #  store actual Pose() object
    #         pdbs.append(pose_from_pdb('benchmark/'+pdb))
    # pose = pdbs[random.randint(0, len(pdbs))] # picke random pose
    native_pose = Pose()
    native_pose.assign(pose)
    resfile = "./design.resfile"
    with open(resfile, "w") as f:
        f.write("ALLAAxc \n")
        f.write("start\n")
    taskf = prs.rosetta.core.pack.task.TaskFactory()
    taskf.push_back(prs.rosetta.core.pack.task.operation.InitializeFromCommandline())
    taskf.push_back(prs.rosetta.core.pack.task.operation.ReadResfile(resfile))
    packer = prs.rosetta.protocols.minimization_packing.PackRotamersMover(scfxn)
    packer.task_factory(taskf)
    taskf.create_task_and_apply_taskoperations(pose)
    packer.apply(pose)
    print("Optimized scfxn score: ", scfxn(pose))
    print("REF15 Score ", ref15(pose))
    # TODO: What defines our loss, for now use REF15 or bloss62 matrix
    bloss62 = substitution_matrices.load("BLOSUM62")
    similar = pairwise2.align.globaldx(
        pose.sequence(), native_pose.sequence(), bloss62, score_only=True
    ) / len(pose.sequence())
    print("similarity ::", similar)

    # return inverse score as we minimize

    return -similar
