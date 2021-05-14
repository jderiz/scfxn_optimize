import logging
import time

import pyrosetta as prs
from pyrosetta import get_fa_scorefxn
import benchmark_prot_fetcher

prs.init(
    options="-ex1 -ex2", set_logging_handler=True, extra_options="-linmem_ig 10 -archive_on_disk /tmp/rosetta"
)  # no output from the design process
prs.logging_support.set_logging_sink()
logger = logging.getLogger("rosetta")
logger.setLevel(logging.ERROR)



global crystals, pre_relaxed

def init():

    # get benchmark protein crystal structures.
    crystals = benchmark_prot_fetcher.get_crystals()
    pre_relaxed = benchmark_prot_fetcher.get_all() 

def relax_with_config(fa_reps):

    # relax the structure and compare to default fast relax
    avg_over = 10
    scfxn  = get_fa_scorefxn()
    prot_name = random.choice(crystals.keys()) 
    pose = crystals[prot_name]
    default_pose = pre_relaxed[prot_name]
    # empty file line vector
    svec = prs.rosetta.std.vector_std_string()

    relax_protocol = prs.rosetta.protocols.relax.FastRelax()
    svec.append("repeat 5")
    svec.append("coord_cst_weight 1.0")
    svec.append("scale:fa_rep "+ fa_reps[0])
    svec.append("repack")
    svec.append("scale:fa_rep " + fa_reps[1])
    svec.append("min 0.01")
    svec.append("coord_cst_weight 0.5")
    svec.append("scale:fa_rep " + fa_reps[2])
    svec.append("repack")
    svec.append("scale:fa_rep " + fa_reps[3]) 
    svec.append("min 0.01")
    svec.append("coord_cst_weight 0.0")
    svec.append("scale:fa_rep " + fa_reps[4])
    svec.append("repack")
    svec.append("scale:fa_rep " + fa_reps[5])
    svec.append("min 0.01")
    svec.append("coord_cst_weight 0.0")
    svec.append("scale:fa_rep " + fa_reps[6])
    svec.append("repack")
    svec.append("min 0.00001")
    svec.append("accept_to_best")
    svec.append("endrepeat")
    
    relax_protocol.set_script_from_lines(svec)
     
    scores = [scfxn(relax_protocol.apply(pose, scfxn))  for _ in range(avg_over)]
    score = sum(scores) / len(scores)
    default_score = scfxn(default_pose)
    return default_score - score

     

