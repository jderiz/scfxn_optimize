#!/usr/bin/env python3
import logging
import os
import random

import pyrosetta as prs
from Bio import pairwise2
from Bio.Align import substitution_matrices
from pyrosetta import (Pose, PyJobDistributor, get_fa_scorefxn, init,
                       pose_from_pdb)
from skopt.utils import use_named_args
import time
import create_scfxn
from hyperparams import scfxn_ref15_space
from setup import runs_per_config


def initialize():
    """initialize pyrosetta
    """
    global prs
    print('INITIALIZED')
    prs.init(
    options="-ex1 -ex2", set_logging_handler=True, extra_options="-linmem_ig 10 -archive_on_disk /tmp/rosetta -mute core -mute"
    )  # no output from the design process
    prs.logging_support.set_logging_sink()
    logger = logging.getLogger("rosetta")
    logger.setLevel(logging.ERROR)
    logger.handlers.pop() # remove stream handler

    rosetta_log_handler = logging.FileHandler('rosetta.log')
    rosetta_log_handler.setLevel(logging.DEBUG)
    logger.addHandler(rosetta_log_handler)
    global pdbs 
    pdbs = {}

    # get all protein pdbs
    for pdb in os.listdir("benchmark"):
        if pdb.endswith(".pdb"):
            #  store actual Pose() object
            pdbs.update({pdb.split('.')[0]: pose_from_pdb("benchmark/" + pdb)})




@use_named_args(dimensions=scfxn_ref15_space)
def design_with_config(**config):
    start_time = time.time()  #  Runtime measuring
    print('DESIGNING')
    ref15 = (
        get_fa_scorefxn()
    )  # REF15 score function with default weights for loss calulation
    scfxn = create_scfxn.creat_scfxn_from_config(
        config=config
    )  # optimization score Function
    
    # pick random
    # pose = random.choice(pdbs)
    # key = random.choice(pdbs.keys())
    # pose = pdbs[key]
    pose = pdbs['1K9P']

    # copy pose for comparison after design
    native_pose = Pose()
    native_pose.assign(pose)
    resfile = "./design.resfile"
    with open(resfile, "w") as f:
        f.write("ALLAAxc \n")
        f.write("start\n")

    # def run(pose):
    taskf = prs.rosetta.core.pack.task.TaskFactory()
    taskf.push_back(prs.rosetta.core.pack.task.operation.InitializeFromCommandline())
    taskf.push_back(prs.rosetta.core.pack.task.operation.ReadResfile(resfile))
    packer = prs.rosetta.protocols.minimization_packing.PackRotamersMover(scfxn)
    packer.task_factory(taskf)
    taskf.create_task_and_apply_taskoperations(pose)
    packer.apply(pose)
    # TODO: What defines our loss, for now use REF15 or bloss62 matrix
    bloss62 = substitution_matrices.load("BLOSUM62")
    # compute normalizes similarity
    similar = pairwise2.align.globaldx(
        pose.sequence(), native_pose.sequence(), bloss62, score_only=True
    ) / len(pose.sequence())

    # moritz says its okay to return energy normalized by length
    result = {"sequence": pose.sequence(), "prot_len": len(pose.sequence()), "bloss62": -similar, "ref15": (ref15(pose)/len(pose.sequence())), "scfxn": (scfxn(pose)/len(pose.sequence()))}
    print('DESIGN_DONE: ',result)
    took = time.time() - start_time
    print("Took: {} to run".format(time.strftime("%H: %M: %S", time.gmtime(took))))

    # This has to be serializable in order to get pickled and send back to parent
    return result
