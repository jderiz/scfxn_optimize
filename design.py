#!/usr/bin/env python3
import random

import pyrosetta as prs
from pyrosetta import *
import create_scfxn
from skopt.utils import use_named_args
from hyperparams import scfxn_ref15_space

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
    ref15 = get_fa_scorefxn()
    scfxn = create_scfxn.creat_scfxn_from_config(config=config)
    pose = pose_from_pdb('benchmark/1K9P_A_relax_0001.pdb')  # easiest struct for optimizer picking
    starting_pose = Pose()
    resfile = "./design.resfile"
    with open(resfile, "w") as f:
        f.write("ALLAAxc \n")
        f.write("start\n")

    taskf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    taskf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    taskf.push_back(pyrosetta.rosetta.core.pack.task.operation.ReadResfile(resfile))
    # taskf.push_back(operation.RestrictToRepacking())
    packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scfxn)
    packer.task_factory(taskf)
    # print(taskf.create_task_and_apply_taskoperations(pose))
    packer.apply(pose)

    print("Optimized scfxn score: ", scfxn(pose))
    print('REF15 Score ', ref15(pose))
    # we optimize the Designing function but eval on standard ref15
    return ref15(pose)


class Designer(object):
    def __init__(self, filename=None, scorefxn=None):
        """Instatiate Designer"""
        if filename:
            self.pose = pose_from_pdb(filename=filename)
        self.starting_pose = Pose()
        self.starting_pose.assign(self.pose)
        self.scorefxn = get_fa_scorefxn()

        if scorefxn:
            self.scorefxn_G = scorefxn
            print("running with custom scorefxn_G")
        else:
            self.scorefxn_G = get_fa_scorefxn()
            print("running with REF15 scorefxn_G")

        # self.relax = prs.rosetta.protocols.relax.FastRelax()

    def fixbb(self, scfxn=None):

        """
            Performs the RosettaDesign protocol to change a structure's
            amino acid sequence while maintaining a FIXED backbone.
            Generates the structure.pdb file.
            """
        if scfxn:
            self.scorefxn_G = scfxn

        resfile = "./design.resfile"
        with open(resfile, "w") as f:
            f.write("ALLAAxc \n")
            f.write("start\n")
        # resfile = prs.rosetta.core.pack.task.operation.ReadResfile("./design.resfile")
        # task = prs.rosetta.core.pack.task.TaskFactory()
        # task.push_back(resfile)

        # job = PyJobDistributor("fixbb", 20, self.scorefxn)
        # job.native_pose = self.starting_pose
        #
        # while not job.job_complete:
        #     self.pose.assign(self.starting_pose)
        #     fixbb.apply(self.pose)
        #     #self.relax.apply(self.pose)
        #     fixbb.output_decoy(self.pose)

        score_fxn = self.scorefxn_G
        pose = self.pose
        taskf = pyrosetta.rosetta.core.pack.task.TaskFactory()
        taskf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
        taskf.push_back(pyrosetta.rosetta.core.pack.task.operation.ReadResfile(resfile))
        # taskf.push_back(operation.RestrictToRepacking())
        packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(score_fxn)
        packer.task_factory(taskf)
        print(taskf.create_task_and_apply_taskoperations(pose))
        packer.apply(pose)

        print(self.scorefxn_G(pose))
        return self.scorefxn_G(pose)

    @use_named_args(dimensions=scfxn_ref15_space)
    @staticmethod
    def test(**config):
        print(config)
        print(create_scfxn.creat_scfxn_from_config(config=config))
        return random.randint(0, 100)


if __name__ == '__main__':
    prs.init("-ex1 -ex2")
    designer = Designer('/home/iwe7/meiler/scfxn_optimize/benchmark/1K9P_A_relax_0001.pdb')
    designer.fixbb()
