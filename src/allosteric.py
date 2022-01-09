import logging
import os
import random
import time

import pandas as pd
import pyrosetta as prs
import ray
from pyrosetta import Pose, get_fa_scorefxn
from pyrosetta.distributed.packed_pose.core import PackedPose


def initialize():

    prs.init(
        options="-ex1 -ex2",
        set_logging_handler=True,
        extra_options="-linmem_ig 10 -archive_on_disk /tmp/rosetta -mute core -mute basic -mute protocols"
    )  # no output from the design process
    prs.logging_support.set_logging_sink()
    logger = logging.getLogger("rosetta")
    logger.setLevel(logging.ERROR)
    # get benchmark protein crystal structures.
    global names
    names = []

    for pdb in os.listdir("../benchmark/allosteric"):
        if pdb.endswith(".pdb"):
            #  store actual Pose() object
            names.append(pdb)


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


def calc_mean_distance(s1, s2):
    summ = 0

    for idx in range(len(s1)):
        summ += calc_distance(s1[idx], s2[idx])

    return summ/len(s1)


def relax_with_config(fa_reps, run, pdb, target=None):
    st = time.time()
    scfxn = get_fa_scorefxn()

    # add ending for correct path
    # pdb_path = pdb+'.pdb'

    ub_dict = {
        "3s0bA.pdb": "3S0A.clean.pdb",
        "1k9kA.pdb": "1K9P.clean.pdb",
        "1f4vA.pdb": "3CHY.clean.pdb",
        "3zjaA.pdb": "3ZK0.clean.pdb",
        "6q21A.pdb": "4q21A.pdb",
        "1avsA.pdb": "1TOP.clean.pdb",
        "1lfaA.pdb": "1MQ9.clean.pdb",
        "1d5wA.pdb": "1d5bA.pdb"
    }

    if target is None:
        unbound: Pose = prs.pose_from_pdb(
            "../benchmark/allosteric/"+ub_dict[pdb])
    elif target:
        unbound: Pose = prs.pose_from_pdb('../benchmark/allosteric/'+target)
    pose: Pose = prs.pose_from_pdb("../benchmark/allosteric/"+pdb)
    # get normalization constants
    start_phi, start_psi = get_phi_psi(pose)
    unb_phi, unb_psi = get_phi_psi(unbound)
    # get Initial Phi/Psi distribution
    phi_norm_const = calc_mean_distance(start_phi, unb_phi)
    psi_norm_const = calc_mean_distance(start_psi, unb_psi)
    start_rmsd = prs.rosetta.core.scoring.CA_rmsd(
        pose, unbound)  # start=start, end=end)  # aligns automatically
    start_ref15 = scfxn(pose)
    # empty file line vector
    svec = prs.rosetta.std.vector_std_string()
    # init relax with score func
    relax_protocol = prs.rosetta.protocols.relax.FastRelax(scfxn)
    # write relax script to vector_std_string

    if fa_reps:  # if config is not Falsey make script with config
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
    relax_protocol.apply(pose)

    rmsd = prs.rosetta.core.scoring.CA_rmsd(
        pose, unbound)  # start=start, end=end)  # aligns automatically

    # TORSION ANGLES
    pose_phi, pose_psi = get_phi_psi(pose)
    phi_deviation = calc_mean_distance(pose_phi, unb_phi)
    psi_deviation = calc_mean_distance(pose_psi, unb_psi)

    # SCORING
    ref15 = scfxn(pose)
    ref15 = ref15/len(pose)
    took = time.strftime("%H:%M:%S", time.gmtime(time.time()-st))
    score = phi_deviation/phi_norm_const \
        + psi_deviation/psi_norm_const \
        + rmsd/start_rmsd
    res = {
        "run": run,
        "pose_phi": pose_phi,
        "pose_psi": pose_psi,
        "target_phi": unb_phi,
        "target_psi": unb_psi,
        "phi_deviation": phi_deviation,
        "psi_deviation": psi_deviation,
        "torsions": phi_deviation+psi_deviation,
        "rmsd": rmsd,
        "config": fa_reps,
        "ref15": ref15,
        "score": score,
        "prot": pdb.split('.')[0],
        "pose": PackedPose(pose),
        "took": took
    }

    return res
