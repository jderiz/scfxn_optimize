import logging
import os
import random
import time

import pyrosetta as prs
from pyrosetta import get_fa_scorefxn
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

    for pdb in os.listdir("benchmark"):
        if pdb.endswith(".pdb"):
            #  store actual Pose() object
            names.append(pdb)


def relax_with_config(pdb, fa_reps):

    # relax the structure and compare to default fast relax
    scfxn = get_fa_scorefxn()
    print('check pdb supplied: ', pdb)
    print('fa_reps: ', fa_reps)

    if pdb:
        # add ending for correct path
        pdb = pdb+'.pdb'
    else:
        pdb = random.choice(names)
    pose = prs.pose_from_pdb("benchmark/crystal/crystal_"+pdb)
    # default_pose = prs.pose_from_pdb("benchmark/1K9P.pdb")
    # empty file line vector
    svec = prs.rosetta.std.vector_std_string()
    # init relax with score func
    relax_protocol = prs.rosetta.protocols.relax.FastRelax(scfxn)
    # write relax script to vector_std_string

    svec.append("repeat 10")
    svec.append("coord_cst_weight 1.0")
    svec.append("scale:fa_rep 0.01")  # + str(fa_reps[0]))
    svec.append("repack")
    # svec.append("scale:fa_rep " + str(fa_reps[1]))
    svec.append("min 0.01")
    svec.append("coord_cst_weight 0.5")
    svec.append("scale:fa_rep 0.33")  # + str(fa_reps[2]))
    svec.append("repack")
    # svec.append("scale:fa_rep " + str(fa_reps[3]))
    svec.append("min 0.01")
    svec.append("coord_cst_weight 0.0")
    svec.append("scale:fa_rep 0.66")  # + str(fa_reps[4]))
    svec.append("repack")
    # svec.append("scale:fa_rep " + str(fa_reps[5]))
    svec.append("min 0.01")
    svec.append("coord_cst_weight 0.0")
    svec.append("scale:fa_rep 0.99")  # + str(fa_reps[6]))
    svec.append("repack")
    svec.append("min 0.00001")
    svec.append("accept_to_best")
    svec.append("endrepeat")
    # make relax use the script
    relax_protocol.set_script_from_lines(svec)
    # evaluate
    print('RUN RELAX')
    relax_protocol.apply(pose)
    score = scfxn(pose)
    # default_score = scfxn(default_pose)
    # normalize by lenght

    res = {
        "score": score/len(pose),
        "prot": pdb.split('.')[0],
        "pose": PackedPose(pose)
    }
    # print(res)

    return res
