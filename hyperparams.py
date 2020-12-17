from skopt.space import Categorical, Real

scfxn_ref15_dimensions = {
    "fa_atr": Real(high=2, low=0),
    'fa_rep': Real(high=2, low=0),
    'fa_sol': Real(high=2, low=0),
    'fa_intra_rep': Real(high=2, low=0),
    'fa_intra_sol_xover4': Real(high=2, low=0),
    'lk_ball_wtd': Real(high=2, low=0),
    'fa_elec': Real(high=2, low=0),
    'pro_close': Real(high=2, low=0),
    'hbond_sr_bb': Real(high=2, low=0),
    'hbond_lr_bb': Real(high=2, low=0),
    'hbond_bb_sc': Real(high=2, low=0),
    'hbond_sc': Real(high=2, low=0),
    'dslf_fa13': Real(high=2, low=0),
    'omega': Real(high=2, low=0),
    'fa_dun': Real(high=2, low=0),
    'p_aa_pp': Real(high=2, low=0),
    'yhh_planarity': Real(high=2, low=0),
    'ref': Real(high=2, low=0),
    'rama_prepro': Real(high=2, low=0),
}
ref15_weights = [('fa_atr', 1),
                 ('fa_rep', 0.55),
                 ('fa_sol', 1),
                 ('fa_intra_rep', 0.005),
                 ('fa_intra_sol_xover4', 1),
                 ('lk_ball_wtd', 1),
                 ('fa_elec', 1),
                 ('pro_close', 1.25),
                 ('hbond_sr_bb', 1),
                 ('hbond_lr_bb', 1),
                 ('hbond_bb_sc', 1),
                 ('hbond_sc', 1),
                 ('dslf_fa13', 1.25),
                 ('omega', 0.4),
                 ('fa_dun', 0.7),
                 ('p_aa_pp', 0.6),
                 ('yhh_planarity', 0.625),
                 ('ref', 1),
                 ('rama_prepro', 0.45)]

relax_init_fa_reps = [0.040, 0.051, 0.265, 0.280, 0.559, 0.581, 1]
relax_init_coord_cst_weight = [1.0, 0.5, 0.0, 0.0]
global _range
_range = 0.25
def get_ref15() -> list:
    return ref15_weights

def get_dimensions() -> list:
    return [Real(high=(value + value*_range),
                              low=(value - value*_range),
                              name=name) for name, value in ref15_weights]


def set_range(new_range: float) -> None:
    global _range

    if new_range < 0 or new_range > 1:
        raise ValueError('range has to be within [0, 1]')
    else:
        _range = new_range
