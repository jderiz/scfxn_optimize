#!/usr/bin/env python3


import logging

from skopt import Optimizer, Space, callbacks


class Opt:
    """
        This class holds the optimizer instance from skopt and manages incorporating new results as         well as supplying the manager with new configurations to test 
    """

    def __init__(self, estimator, dimensions, cooldown=False, **kwargs):
        self.opti: Optimizer = Optimizer(
            base_estimator=estimator, dimensions=dimensions, **kwargs)
        self.logger.log("initialized with ")
        self.cooldown = cooldown
        self.n_calls = 0

        # LOOKUP TABLE FOR COOLDOWN
        self.xi_kappa_lookup = pd.DataFrame(None, index=range(n_calls))
        self.xi_kappa_lookup["iter"] = range(1, n_calls + 1)
        self.xi_kappa_lookup["geospace"] = np.geomspace(0.001, 1, num=n_calls)
        self.xi_kappa_lookup["xi"] = xi - \
            xi_kappa_lookup.geospace * (xi - final_xi)
        self.xi_kappa_lookup["kappa"] = kappa - \
            xi_kappa_lookup.geospace * (kappa - final_kappa)

    def handle_result(config, res) -> bool:
        self.opti.tell(config, res)

        if self.cooldown:
            sefl.cooldown()

        return True

    def get_next_config() -> list:

        return optimizer.ask()

    def cooldown() -> None:

        opti
        pass
