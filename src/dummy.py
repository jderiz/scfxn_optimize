import logging
import random
import time

import ray

logger = logging.getLogger('Dummy')
logger.setLevel(logging.DEBUG)

# @ray.remote


def dummy_objective(config, run, pdb) -> dict:
    time.sleep(random.randint(2, 3))

    return {
        "bloss62": random.randint(1, 100),
        "ref15": random.randint(1, 50),
        "scfxn": random.randint(1, 46),
        "score": random.randint(1, 20),
        "rmsd": random.randint(1, 5),
        "run": run,
        "pdb": pdb,
        "config": config

    }
