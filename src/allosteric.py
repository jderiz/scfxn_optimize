import logging
import math
import os
import random
import time

import pandas as pd
import pyrosetta as prs
import ray
from pyrosetta import Pose, get_fa_scorefxn
from pyrosetta.distributed.packed_pose.core import PackedPose

logger = logging.getLogger("rosetta")
logger.setLevel(logging.DEBUG)


def initialize():

    prs.init(options="-ex1 -ex2",
             set_logging_handler=True,
             extra_options="-linmem_ig 10 -archive_on_disk /tmp/rosetta -mute core -mute basic -mute protocols"
             )  # no output from the design process
    prs.logging_support.set_logging_sink()


def get_phi_psi(pose):
    phi_angles = []
    psi_angles = []

    for residx in range(1, pose.size()+1):
        phi = pose.phi(residx)
        psi = pose.psi(residx)
        phi_angles.append(phi)
        psi_angles.append(psi)

    return(phi_angles, psi_angles)


def calc_distance(alpha, beta):
    phi = abs(beta-alpha) % 360

    if phi > 180:
        return 360-phi
    else:
        return phi

# see @Sampling the conformational space of allosteric transition using LoopHash
# For Torsion RMSD definition


def calc_torsion_rmsd(s1, s2):
    # if not equal lenght then stop, we dont know where missing residues is
    # TODO: implement proper error handling later
    s1phi, s1psi = get_phi_psi(s1)
    s2phi, s2psi = get_phi_psi(s2)

    resnums = len(s1)
    phi_summ = 0
    psi_summ = 0

    for idx in range(resnums):
        phi_deviation = calc_distance(s1phi[idx], s2phi[idx])
        psi_deviation = calc_distance(s1psi[idx], s2psi[idx])
        phi_summ += phi_deviation**2
        psi_summ += psi_deviation**2

    return math.sqrt((phi_summ+psi_summ)/(2*resnums))


def relax_with_config(fa_reps, run, pdb, target):
    # get start time for timing
    st = time.time()
    # load REF15 scorefunction
    scfxn = get_fa_scorefxn()
    os.chdir("/home/iwe7/scfxn_optimize/src/")
    # dictionary storing the benchmark pairs TODO: deprecated
    ub_dict = {
        "3s0bA.pdb": "3S0A.pdb",
        "1k9kA.pdb": "1K9P.pdb",
        "1f4vA.pdb": "3CHY.pdb",
        "3zjaA.pdb": "3ZK0.pdb",
        "6q21A.pdb": "4Q21.pdb",
        "1lfaA.pdb": "1MQ9.pdb",
        "1d5wA.pdb": "1D5B.pdb"
    }
    logger.debug('IN: %s', os.getcwd())
    target_pose: Pose = prs.pose_from_pdb(target)
    work_pose: Pose = prs.pose_from_pdb(pdb)

    # get Initial rmsd of Phi/Psi angles between the bound and unbound state
    torsion_norm_const = calc_torsion_rmsd(work_pose, target_pose)

    start_rmsd = prs.rosetta.core.scoring.CA_rmsd(
        work_pose, target_pose)  # aligns automatically and return C-alpha RMSD
    start_ref15 = scfxn(work_pose)  # get initial REF15 score for reference
    # empty file line vector
    svec = prs.rosetta.std.vector_std_string()
    # init relax with score func
    relax_protocol = prs.rosetta.protocols.relax.FastRelax(scfxn)
    # write relax script to vector_std_string

    if fa_reps:  # if config supplied make script with config
        svec.append("repeat 5")
        svec.append("coord_cst_weight 1.0")
        svec.append("scale:fa_rep " + str(fa_reps[0]))
        svec.append("repack")
        svec.append("scale:fa_rep " + str(fa_reps[1]))
        svec.append("min 0.01")
        svec.append("coord_cst_weight 0.5")
        svec.append("scale:fa_rep " + str(fa_reps[2]))
        svec.append("repack")
        svec.append("scale:fa_rep " + str(fa_reps[3]))
        svec.append("min 0.01")
        svec.append("coord_cst_weight 0.0")
        svec.append("scale:fa_rep " + str(fa_reps[4]))
        svec.append("repack")
        svec.append("scale:fa_rep " + str(fa_reps[5]))
        svec.append("min 0.01")
        svec.append("coord_cst_weight 0.0")
        svec.append("scale:fa_rep " + str(fa_reps[6]))
        svec.append("repack")
        svec.append("min 0.00001")
        svec.append("accept_to_best")
        svec.append("endrepeat")
        # make relax use the script
        relax_protocol.set_script_from_lines(svec)
    # RUN RELAX
    # relax_protocol.apply(work_pose)

    rmsd = prs.rosetta.core.scoring.CA_rmsd(
        work_pose, target_pose)  # start=start, end=end)  # aligns automatically

    # TORSION RMSD AFTER RELAX
    torsion_rmsd = calc_torsion_rmsd(work_pose, target_pose)
    # SCORING
    ref15 = scfxn(work_pose)
    ref15 = ref15/len(work_pose.residues)
    took = time.strftime("%H:%M:%S", time.gmtime(time.time()-st))
    score = (torsion_rmsd/torsion_norm_const + rmsd/start_rmsd)/2
    res = {
        "run": run,
        "torsion_rmsd": torsion_rmsd,
        "rmsd": rmsd,
        "config": fa_reps,
        "ref15": ref15,
        "score": score,
        "pose": PackedPose(work_pose),
        "took": took
    }

    return res
