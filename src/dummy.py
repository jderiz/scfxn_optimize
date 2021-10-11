import random
import time
import logging

logger = logging.getLogger('Dummy')

def dummy_objective(config, run, pdb) -> dict:
    time.sleep(random.randint(2,3))
    logger.debug('EVAL: %s:::: RUN: %s', config, run)
    return {
        "bloss62": random.randint(1, 100),
        "ref15": random.randint(1, 50),
        "scfxn": random.randint(1, 46),
        "score": random.randint(1, 20),
        "rmsd": random.randint(1,5),
        "run": run, 
        "pdb": pdb, 
        "config": config
        
    }
