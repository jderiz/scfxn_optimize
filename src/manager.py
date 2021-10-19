import logging
import os
import pickle
import sys
import threading
import time
from functools import partial

import numpy as np
import pandas as pd
import ray

import config
from bayesopt import BayesOpt
from parallel import Distributor

logger = logging.getLogger('OptimizationManager')
logger.setLevel(logging.DEBUG)
lock = threading.Lock()


class Singleton(type):
    """docstring for Singleton"""
    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            with lock:
                if cls not in cls._instances:
                    cls._instances[cls] = super(
                        Singleton, cls).__call__(*args, **kwargs)

        return cls._instances[cls]


@ray.remote
class OptimizationManager():

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
        self.identify = identifier
        self.base_estimator = estimator
        self.evals = evals
        self._DONE = False
        self.pandas = save_pandas
        self.loss = loss
        self.evals_done = 0
        self.pdb = pdb
        self.n_cores = n_cores
        self.rpc = rpc
        self.test_run = test_run
        # TEST CASE

        if test_run:
            logger.info('TEST RUN')
            self.objective = config._dummy_objective
            self.loss_value = 'ref15'
            self.init_method = None
            self.pandas = False
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
            manager_callback=self.log_res_and_update, hpc=hpc, cpus=n_cores, initializer=self.init_method)

        # BOOKKEEPING
        self.results = pd.DataFrame()
        logger.debug('initialized OptimizationManager')

    def log_res_and_update(self, map_res) -> None:
        # TODO: find type of future.result() in dask
        self.evals_done += 1
        # logger.debug('%s', map_res)
        self.results = self.results.append(
            pd.DataFrame(map_res), ignore_index=True)
        self.optimizer.update_prior(
            map_res[0]['config'], sum(x[self.loss_value] for x in map_res)/len(map_res))

        if self.evals_done < self.evals:
            self.make_batch()
        elif self.evals_done >= self.evals:
            self._save_and_exit()

    def make_batch(self):
        config = self.optimizer.get_next_config()
        self.distributor.distribute(func=self.objective,
                                    params=config,
                                    pdb=self.pdb,
                                    num_workers=self.rpc,
                                    run=self.evals_done+1)

    def _save_and_exit(self) -> bool:
        """
            Saves the results stored in the DataFrame and reports to console
        """
        logger.debug(' ')

        if not self.test_run:
            if self.pandas:
                # save pandas DataFrame with correct column names if already DataFrame then nothing changes
                df = pd.DataFrame(self.results)
                weights = df.config.apply(lambda x: pd.Series(x))
                weights.columns = ["fa_rep_" + str(num) for num in range(7)]
                results = pd.concat([df, weights], axis=1)
            with open(
                "results/{}_{}_res_{}.pkl".format(self.identify,
                                                  self.base_estimator, self.evals),
                "wb",
            ) as file:
                pickle.dump(self.results, file)
        else:
            print(self.results.head())
        self._DONE = True
        self.distributor.terminate()
        print("TERMINATING")

    def run(self) -> None:
        """
            Runs the Manager and 
        """
        logger.debug('RUN OptimizationManager')
        # map initial runs workers/rpc rpc times

        for _ in range(int(self.n_cores/self.rpc)):
            self.distributor.distribute(func=self.objective,
                                        params=self.optimizer.get_next_config(),
                                        pdb=self.pdb,
                                        num_workers=self.rpc,
                                        run=self.evals_done+1)

        # logger.debug("distributor.jobs. len %s",
        #              self.distributor.jobs.__len__())
        # logger.debug("job.wait() result %s",
        #              self.distributor.jobs.pop().wait())

    def is_done(self):
        return self._DONE


_manager = OptimizationManager.options(name='manager').remote()
