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
    )  # no output from rosetta
    prs.logging_support.set_logging_sink()
    logger = logging.getLogger("rosetta")
    logger.setLevel(logging.ERROR)
    logger.info('SETUP')
    # get benchmark protein crystal structures.
    global names
    names = []

    # for pdb in os.listdir("benchmark"):
    #     if pdb.endswith(".pdb"):
    #         #  store actual Pose() object
    #         names.append(pdb)
    global wild_type_pose
    print('load WT')
    wild_type_pose = prs.pose_from_pdb("benchmark/RMSD/2lzm.pdb")

    # write relax script to vector_std_string
    print('setup r_select')
    global residue_selector
    residue_selector = prs.rosetta.core.select.residue_selector.ResidueIndexSelector() 
    residue_selector.set_index_range(102, 126)
    
    # make movemap for relax_protocol
    print('setup move_map')
    global move_map
    move_map = prs.MoveMap()
    move_map.set_bb_true_range(102, 126)
    move_map.set_chi(True)
    
    # align based on everything except modelled section
    print('setup algin select')
    global align_selector
    align_selector = prs.rosetta.core.select.residue_selector.NotResidueSelector()
    align_selector.set_residue_selector(residue_selector)
    
    # rmsd metric for scoring, only scores the modelled subset
    print('setup rmsd')
    global rmsd_metric
    rmsd_metric = prs.rosetta.core.simple_metrics.metrics.RMSDMetric()
    rmsd_metric.set_residue_selector(residue_selector)
    rmsd_metric.set_comparison_pose(wild_type_pose)

    # fill position vector for superimpose
    print('make si_vec')
    global si_vec
    si_vec = prs.rosetta.utility.vector1_unsigned_long()
    for i in range(1, 102):
        si_vec.append(i)
    for i in range(127, 165):
        si_vec.append(i)
def relax_with_config(pdb, fa_reps):

    # relax the structure and compare to default fast relax
    scfxn = get_fa_scorefxn()
    # print('check pdb supplied: ', pdb)
    print('fa_reps: ', fa_reps)

    # if pdb:
    #     # add ending for correct path
    #     pdb = pdb+'.pdb'
    # else:
    #     pdb = random.choice(names)
    mutate_pose = prs.pose_from_pdb("benchmark/RMSD/2lc9.pdb")
    mutate_pose_raw = prs.pose_from_pdb("benchmark/RMSD/2lc9.pdb")

    # init relax with score func
    relax_protocol = prs.rosetta.protocols.relax.FastRelax(scfxn)
    global rmsd_metric
    global si_vec
    global move_map
    global wild_type_pose
    # empty file line vector
    svec = prs.rosetta.std.vector_std_string()
    svec.append("repeat 5")
    svec.append("coord_cst_weight 1.0")
    svec.append("scale:fa_rep "+ str(fa_reps[0]))
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
    svec.append("scale:fa_rep "+ str(fa_reps[6]))
    svec.append("repack")
    svec.append("min 0.00001")
    svec.append("accept_to_best")
    svec.append("endrepeat")
    # make relax use the script
    relax_protocol.set_script_from_lines(svec)
    relax_protocol.set_movemap(move_map)
    # print(relax_protocol.get_movemap().show())
   
    # evaluate
    print('RUN RELAX')
    relax_protocol.apply(mutate_pose)

    print('ALIGN STRUCT')
    # align poses based on every residue except the relaxed section
    prs.rosetta.protocols.toolbox.pose_manipulation.superimpose_pose_on_subset_CA(mutate_pose, wild_type_pose, si_vec)

    print('SCORE')
    score = scfxn(mutate_pose)
    rmsd_metric.apply(mutate_pose)
    # default_score = scfxn(default_pose)
    # normalize by lenght

    res = {
        "score": score/len(mutate_pose.residues),
        # "prot": pdb.split('.')[0],
        "pose": PackedPose(mutate_pose),
        "rmsd": mutate_pose.scores['rmsd']
    }
    print(res)

    return res
