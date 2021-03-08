import logging
import time

import pyrosetta as prs
import benchmark_prot_fetcher
prs.init(
    options="-ex1 -ex2", set_logging_handler=True, extra_options="-linmem_ig 10 -archive_on_disk /tmp/rosetta"
)  # no output from the design process
prs.logging_support.set_logging_sink()
logger = logging.getLogger("rosetta")
logger.setLevel(logging.ERROR)


def relax_with_config(fa_reps):
    script_manager = prs.rosetta.protocols.relax.RelaxScriptManager()
    

    # get benchmark protein
    pdbs = benchmark_prot_fetcher.get_crystal()


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
