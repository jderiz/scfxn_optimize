import logging

import numpy as np
import pandas as pd
from skopt import Optimizer, Space, callbacks
import config   

class BayesOpt:
    """
        This class holds the optimizer instance and manages incorporating new results 
        as well as supplying the manager with new configurations to test 
    """

    def __init__(self, estimator='RF', **kwargs):
        self.logger = logging.getLogger('BayesOpt')
        self.logger.setLevel(logging.DEBUG)

        # 
        try:
            dimensions = Space.from_yaml(kwargs['dimensions'])
        except FileNotFoundError as e:
            print(e)
            print(
                "No space file has been defined the optimizer cannot be defined"
            )
            exit()
        self.opti: Optimizer = Optimizer(
            base_estimator=estimator, dimensions=dimensions)
        self.logger.info("Optimizer initialized")
        self.cooldown = kwargs['cooldown']
        self.evals = kwargs['evals']
        self.calls = 0

        # LOOKUP TABLE FOR COOLDOWN
        self.xi_kappa_lookup = pd.DataFrame(None, index=range(self.evals))
        self.xi_kappa_lookup["iter"] = range(1, self.evals + 1)
        self.xi_kappa_lookup["geospace"] = np.geomspace(
            0.001, 1, num=self.evals)
        self.xi_kappa_lookup["xi"] = config.xi - \
            self.xi_kappa_lookup.geospace * (config.xi - config.final_xi)
        self.xi_kappa_lookup["kappa"] = config.kappa - \
            self.xi_kappa_lookup.geospace * \
            (config.kappa - config.final_kappa)



    def handle_result(self, config, res) -> bool:
        self.opti.tell(config, res)

        if self.cooldown:
            self._cooldown()

        return True

    def get_next_config(self) -> list:
        return self.opti.ask()

    def _cooldown(self) -> None:

        new_kappa = self.xi_kappa_lookup.kappa[self.calls]
        new_xi = self.xi_kappa_lookup.xi[self.calls]
        print("UPDATE XI KAPPA \n", new_xi, new_kappa)
        self.opti.acq_optimizer_kwargs.update(
            {"xi": new_xi, "kappa": new_kappa})
        # update aquisition
        self.opti.update_next()

    def update_prior(self, x, y):
        """
        updates the Optimizers prior 
        """
        self.opti.tell(x, y)
