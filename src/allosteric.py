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
        options="-ex1 -ex2", set_logging_handler=True, extra_options="-linmem_ig 10 -archive_on_disk /tmp/rosetta -mute core -mute basic -mute protocols"
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


def relax_with_config(fa_reps, run, pdb):
    st = time.time()
    scfxn = get_fa_scorefxn()

    # add ending for correct path
    # pdb_path = pdb+'.pdb'

    ub_dict = {
        "3S0B.clean.pdb": "3S0A.clean.pdb",
        "1k9kA.pdb": "1K9P.clean.pdb",
        "1f4vA.pdb": "3CHY.clean.pdb",
        "3zjaA.pdb": "3ZK0.clean.pdb",
        # "6q21A": "4Q21.clean.pdb",
        # "1avsA": "1TOP.clean.pdb",
        # "1lfaA": "1MQ9.clean.pdb",
        # "1d5wA": "1DBW.clean.pdb"
    }
    unbound: Pose = prs.pose_from_pdb(
        "../benchmark/allosteric/"+ub_dict[pdb])
    # print(os.getcwd())
    pose: Pose = prs.pose_from_pdb("../benchmark/allosteric/"+pdb)
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

    # get shorter sequence and use for start - end

    # if (len(pose.residues) >= len(unbound.residues)):
    #     start = unbound.chain_begin(1)
    #     end = unbound.chain_end(1)
    # else:
    #     start = pose.chain_begin(1)
    #     end = pose.chain_end(1)
    rmsd = prs.rosetta.core.scoring.CA_rmsd(
        pose, unbound)  # start=start, end=end)  # aligns automatically

    # TORSION ANGLES
    pose_phi, pose_psi = get_phi_psi(pose)
    unb_phi, unb_psi = get_phi_psi(unbound)
    pp = pd.DataFrame({"phi": pose_phi, "psi": pose_psi})
    upp = pd.DataFrame({"phi": unb_phi, "psi": unb_psi})
    phi_deviation = (pp.phi - upp.phi).abs().mean()
    psi_deviation = (pp.psi - upp.psi).abs().mean()

    # SCORING
    ref15 = scfxn(pose)
    ref15 = ref15/len(pose)
    took = time.strftime("%H:%M:%S", time.gmtime(time.time()-st))
    score = (phi_deviation+psi_deviation)/2 + (1+rmsd)**2 + ref15
    res = {
        "run": run,
        "pose_phi_psi": pp,
        "target_phi_psi": upp,
        "rmsd": rmsd,
        "config": fa_reps,
        "ref15": ref15,
        "score": score,
        "prot": pdb.split('.')[0],
        "pose": PackedPose(pose),
        "took": took
    }

    return res
