import logging
import time

import pyrosetta as prs
from pyrosetta import get_fa_scorefxn
from pyrosetta.rosetta.protocols.relax import FastRelax
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

    relax_protocol = FastRelax()
    with open('relax_script', 'wb') as f:
        f.write("repeat 5 \n\
                coord_cst_weight 1.0\n\
                scale:fa_rep {}\n\
                repack\n\
                scale:fa_rep {}\n\
                min 0.01\n\
                coord_cst_weight 0.5\n\
                scale:fa_rep {}\n\
                repack\n\
                scale:fa_rep {}\n\
                min 0.01\n\
                coord_cst_weight 0.0\n\
                scale:fa_rep {}\n\
                repack\n\
                scale:fa_rep {}\n\
                min 0.01\n\
                coord_cst_weight 0.0\n\
                scale:fa_rep {}\n\
                repack\n\
                min 0.00001\n\
                accept_to_best\n\
                endrepeat".format(*fa_reps))
     
    scores = [scfxn(relax_protocol.apply(pose, scfxn))  for _ in range(avg_over)]
    score = sum(scores) / len(scores)
    default_score = scfxn(default_pose)
    return default_score - score

     

