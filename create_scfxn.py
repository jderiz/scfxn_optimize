from pyrosetta import ScoreFunction, get_fa_scorefxn

from pyrosetta.rosetta.core.scoring import score_type_from_name


def creat_scfxn_from_config(config):
    scfxn = ScoreFunction()
    print(config)
    for n, w in config.items():
        print(w)
        score_type = score_type_from_name(n)
        scfxn.set_weight(score_type, w)

    return scfxn


def create_ref15():
    return get_fa_scorefxn()
