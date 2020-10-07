from skopt.space import Real, Categorical

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
    # 'problem_type': Categorical()
}

scfxn_ref15_space = [Real(high=5, low=0, name=name) for name, _ in scfxn_ref15_dimensions.items()]
