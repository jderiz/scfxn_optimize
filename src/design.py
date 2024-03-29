#!/usr/bin/env python3
import logging
import os
import random
import time

import pyrosetta as prs
from Bio import pairwise2
from Bio.Align import substitution_matrices
from pyrosetta import (Pose, PyJobDistributor, get_fa_scorefxn, init,
                       pose_from_pdb)
from pyrosetta.distributed.packed_pose.core import PackedPose
from pyrosetta.rosetta.utility import vector1_bool
from skopt.utils import use_named_args

import create_scfxn
import hyperparams

pdbs_dict = {}

def initialize(pdbs, pdb_path):
    """initialize pyrosetta
    """
    global prs
    print('INITIALIZED')
    prs.init(
        options="-ex1 -ex2", set_logging_handler=True, extra_options="-linmem_ig 10 -archive_on_disk /tmp/rosetta -mute core -mute basic -mute protocols"
    )  # no utput from the design process
    # Setup Logging
    prs.logging_support.set_logging_sink()
    logger = logging.getLogger("rosetta")
    logger.setLevel(logging.ERROR)
    # logger.handlers.pop()  # remove stream handler

    rosetta_log_handler = logging.FileHandler('rosetta.log')
    rosetta_log_handler.setLevel(logging.DEBUG)
    logger.addHandler(rosetta_log_handler)

    # Setup proteins and their pssm matrices
    # Proteins
    for pdb in pdbs:
        pdbs_dict.update({pdb.split()[0]: prs.pose_from_pdb(pdb_path+'/'+pdb)})
    # PSSMS
    global pssms
    pssms = {}

    for pssm in os.listdir('PSSMS'):
        # for all benchmark pdbs set the native pose as comparison_pose
        # set ReturnResidueSubsetSelector with whole sequence
        names = pssm.split('_')[1]

        for prot_name in pdbs.keys():
            if prot_name in names:
                # print(prot_name, names)
                seq_len = len(pdbs[prot_name].sequence())
                sub_vec = vector1_bool()

                # except 6Q21 --> 168

                if prot_name == '6Q21':
                    for _ in range(168):
                        sub_vec.append(1)
                else:
                    for _ in range(seq_len):
                        # make  as long as sequence
                        sub_vec.append(1)
                srm = prs.rosetta.protocols.analysis.simple_metrics.SequenceRecoveryMetric()
                # get ResidueSelector that selects entire Sequence
                rselector = prs.rosetta.core.select.residue_selector.ReturnResidueSubsetSelector(
                    sub_vec)
                srm.load_pssm(os.getcwd()+'/PSSMS/' + pssm)
                # compare to native pose
                srm.set_comparison_pose(pdbs[prot_name])
                srm.set_residue_selector(rselector)
                # uses same as set_residue_selector
                srm.set_residue_selector_ref(None)
                pssms.update({prot_name: srm})
        # BLOSS TODO: let blossum compute from prs
        # bloss = prs.rosetta.proticils.analysis.simple_metrics.SequenceSimilarityMetric()

    # print(pdbs, pssms)
    # return pssms


def design_with_config(config) -> dict:
    start_time = time.time()  # Runtime measuring
    print('DESIGNING')
    ref15 = get_fa_scorefxn()  # REF15
    if config == 'ref15':
        scfxn = ref15
    else:
        scfxn = create_scfxn.creat_scfxn_from_config(
            config=config
        )  # optimization score Function

    # pick random protein from pdbs_dict
    prot_name = random.choice(pdbs_dict.keys())
    pose = pdbs_dict[prot_name]


    # copy pose for comparison after design
    native_pose = Pose()
    native_pose.assign(pose)
    resfile = "./design.resfile"
    with open(resfile, "w") as f:
        f.write("ALLAAxc \n")
        f.write("start\n")

    taskf = prs.rosetta.core.pack.task.TaskFactory()
    taskf.push_back(
        prs.rosetta.core.pack.task.operation.InitializeFromCommandline())
    taskf.push_back(prs.rosetta.core.pack.task.operation.ReadResfile(resfile))
    packer = prs.rosetta.protocols.minimization_packing.PackRotamersMover(
        scfxn)
    packer.task_factory(taskf)
    taskf.create_task_and_apply_taskoperations(pose)
    packer.apply(pose)
    # TODO: What defines our loss, for now use REF15 or bloss62 matrix
    bloss62 = substitution_matrices.load("BLOSUM62")
    # compute normalizes similarity
    similar = pairwise2.align.globaldx(
        pose.sequence(), native_pose.sequence(), bloss62, score_only=True
    ) / len(pose.sequence())
    # compute Rosetta SimpleMetrics PSSM
    pssm_score = pssms[prot_name].calculate(pose)
    print('scored with pssm ')

    # moritz says its okay to return energy normalized by length
    # check if pose can be pickled fast and returned

    took = time.time() - start_time

    # This has to be serializable in order to get pickled and send back to parent
    result = {"sequence": pose.sequence(),
              # "pose": PackedPose(pose),
              "prot_len": len(pose.sequence()),
              "prot_name": prot_name,
              "bloss62": -similar,
              "ref15": (ref15(pose)/len(pose.sequence())),
              "scfxn": (scfxn(pose)/len(pose.sequence())),
              "pssm": -pssm_score,
              "runtime": took
              }

    print('DESIGN_DONE: ')
    took = time.time() - start_time
    print("Took: {} to run design on length {}".format(
        time.strftime("%H: %M: %S", time.gmtime(took)), len(pose.sequence())))


    return result
