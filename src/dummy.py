import random
def dummy_objective(config, pdb) -> dict:

    return {
        "bloss62": random.randint(1, 100),
        "ref15": random.randint(1, 50),
        "scfxn": random.randint(1, 46),
        "score": random.randint(1, 20),
        "rmsd": random.randint(1,5)
    }
