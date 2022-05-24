import logging
import os

import numpy as np
import pandas as pd
import ray
from skopt import Optimizer, Space, callbacks

import config


class BayesOpt:
    """
        This class holds the optimizer instance and manages incorporating new results 
        as well as supplying the manager with new configurations to test 
    """

    def __init__(self):
        """
        DUMMY INIT
        """
        pass

    def init(self, dimensions=None, cooldown=True,
             evals=12, base_estimator='RF',
             random_state=5,
             acq_func_kwargs=None,
             n_initial_points=20):
        self.logger = logging.getLogger('BayesOpt')
        self.logger.setLevel(logging.DEBUG)
        self.logger.debug('CWD %s', os.getcwd())
        #
        try:
            dimensions = Space.from_yaml(config.space_dimensions)
        except FileNotFoundError as e:
            print(e)
            print(
                "No space file has been defined the optimizer cannot be defined"
            )
            exit()
        self.opti: Optimizer = Optimizer(
            base_estimator=base_estimator,
            n_initial_points=n_initial_points,
            dimensions=dimensions,
            random_state=random_state,
            acq_func="gp_hedge",
            acq_func_kwargs=acq_func_kwargs)
        self.cooldown = cooldown
        self.evals = evals
        self.updates = 0

        # LOOKUP TABLE FOR COOLDOWN
        self.xi_kappa_lookup = pd.DataFrame(None, index=range(1, self.evals+1))
        self.xi_kappa_lookup["iter"] = range(1, self.evals + 1)
        self.xi_kappa_lookup["geospace"] = np.geomspace(
            0.001, 1, num=self.evals)
        self.xi_kappa_lookup["xi"] = config.xi - \
            self.xi_kappa_lookup.geospace * (config.xi - config.final_xi)
        self.xi_kappa_lookup["kappa"] = config.kappa - \
            self.xi_kappa_lookup.geospace * \
            (config.kappa - config.final_kappa)
        self.logger.info('INITIALIZED')
        self.logger.debug('LOOKUP TABLE \n %s', self.xi_kappa_lookup.head())

    def handle_result(self, config, res) -> bool:
        self._update_prior(config, res)

        if self.cooldown:
            self._cooldown()

        return True

    def get_next_config(self) -> list:

        return self.opti.ask()

    def _cooldown(self) -> None:
        self.logger.debug(self.updates)
        new_kappa = self.xi_kappa_lookup.kappa[self.updates]
        new_xi = self.xi_kappa_lookup.xi[self.updates]
        self.logger.info(
            "UPDATE XI KAPPA \n xi: %f \n kappa: %f", new_xi, new_kappa)
        self.opti.acq_optimizer_kwargs.update(
            {"xi": new_xi, "kappa": new_kappa})
        # update aquisition
        self.opti.update_next()

    def _update_prior(self, x, y):
        """
        updates the Optimizers prior 
        """
        self.opti.tell(x, y)
        self.updates += 1
