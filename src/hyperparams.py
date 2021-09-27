from skopt.space import Categorical, Real

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

# relax_dims = [Real(high=0.051, low=0), Real(high=0.265, low=0.040), Real(high=0.280, low=0.051), Real(high=0.559, low=0.265), Real(low=0.280, high=0.581), Real(low=0.559, high=1), Real(low=0.581, high=1)]
relax_dims = [Real(low=0, high=1),Real(low=0, high=1),Real(low=0, high=1),Real(low=0, high=1),Real(low=0, high=1),Real(low=0, high=1),Real(low=0, high=1)]
relax_init_coord_cst_weight = [1.0, 0.5, 0.0, 0.0]
global _range
_range = 0.25
def get_dimensions() -> list:
    return [Real(high=(value + value*_range),
                              low=(value - value*_range),
                              name=name) for name, value in ref15_weights]


