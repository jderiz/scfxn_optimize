import logging
import os
import pickle
import random
import sys
import time
from functools import partial

import numpy as np
import pandas as pd

import config
from bayesopt import BayesOpt
from parallel import Distributor


class OptimizationManager:
    # dummy init
    def __init__(self):
        pass
    # INIT

    def init(self,
             loss,
             pdb=None,
             estimator="RF",  # "dummy" for random search
             identifier=None,  # string to identify optimization run
             test_run=False,  # if test run eval dummy_objective instead of real
             evals=200,  # configuration evaluations on the objective
             rpc=8,  # runs_per_config n_calls/rpc = evals
             out_dir="results",  # where results get saved
             warm_start=None,  # continue previous optimization run
             n_cores=None,
             mtpc=None,  # maxtasksperchild
             cooldown=True,  # cooldown exploration to exploitation
             space_dimensions=None,  # yaml file with optimizer dimensions
             save_pandas=True,
             hpc=False  # wether to use single node multiprocessing or some hpc manager
             ):
        self.logger = logging.getLogger('OptimizationManager')
        self.identify = identifier
        self.base_estimator = estimator
        self.evals = evals
        self._DONE = False
        self.pandas = save_pandas
        self.loss = loss
        self.batches_done = 0
        self.n_batches = evals/rpc
        self.pdb = pdb
        self.n_cores = n_cores
        self.rpc = rpc
        # TEST CASE

        if test_run:
            print("DUMMY OBJECTIVE")
            objective = dummy_objective
            loss_value = 'ref15'
            init_method = None
            pandas = False
        else:
            self.objective = config._objective
            self.init_method = config._init_method

        # SEARCH SPACE

        # OPTIMIZER
        self.optimizer = BayesOpt(
            random_state=5,
            dimensions=config.space_dimensions,
            base_estimator=estimator,
            acq_func_kwargs=config.acq_func_kwargs,
            n_initial_points=n_cores/rpc,
            cooldown=cooldown,
            evals=evals
        )
        # DISTRIBUTOR
        self.distributor = Distributor(
            callback=self.log_res_and_update, hpc=hpc, workers=n_cores)

        # BOOKKEEPING
        self.results = pd.DataFrame()
        self.logger.debug('initialized OptimizationManager')

    def log_res_and_update(self, config, run, map_res: dict) -> None:
        # TODO: find type of future.result() in dask

        print(map_res)

        for res in map_res:
            res.update({"config": config})
            res.update({"run": run})
        self.results.extend(map_res)
        self.optimizer.update_prior(config, (sum([x[self.loss]
                                                  for x in map_res]) / len(map_res)))

        if self.batches_done < self.n_batches:
            self.make_batch()
        elif self.batches_done == self.n_batches:
            self.save_and_exit()

    def run(self) -> None:

        # map initial runs workers/rpc rpc times

        for _ in range(self.n_cores/self.rpc):
            self.distributor.distribute(func=partial(
                self.objective, pdb=self.pdb),
                params=self.optimizer.get_next_config(),
                num_workers=self.rpc,
                run=self.batches_done+1)
        self.batches_done+=1
        # wait in this loop until DONE

        while not self._DONE:
            time.sleep(5)

    def make_batch(self):
        config = self.optimizer.get_next_config()
        self.distributor.distribute(func=partial(
            self.objective, self.pdb),
            params=config,
            num_workers=self.rpc,
            run=self.batches_done+1)
        self.batches_done+=1
    def save_and_exit(self) -> bool:

        if self.pandas:
            # save pandas DataFrame with correct column names
            df = pd.DataFrame(results)
            weights = df.config.apply(lambda x: pd.Series(x))
            weights.columns = ["fa_rep_" + str(num) for num in range(7)]
            results = pd.concat([df, weights], axis=1)
        with open(
            "results/{}_{}_res_{}.pkl".format(self.identify,
                                              self.base_estimator, self.evals),
            "wb",
        ) as file:
            pickle.dump(results, file)
        print("TERMINATING")
        self._DONE = True

    @ staticmethod
    def no_optimize(identify, config_path, pdb, evals):
        pass

# use for test purposes


def dummy_objective(pdb, config) -> dict:
    print('TEST RUN')
    # time.sleep(random.randint(5, 15))

    return {
        "bloss62": random.randint(1, 100),
        "ref15": random.randint(1, 50),
        "scfxn": random.randint(1, 46),
        "score": random.randint(1, 20),
    }


# SINGLETON
Manager = OptimizationManager()
