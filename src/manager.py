import asyncio
import logging
import os
import pickle
import sys
import threading
import time

import numpy as np
import pandas as pd
import ray
from ray.util import inspect_serializability

import config

# from bayesopt import BayesOpt
# from parallel import Distributor


@ray.remote
class OptimizationManager():
    """
    Optimization Manager is a single Actor that handles the communication between optimizer and distributor actor and does all the saving to file
    """
    # dummy init

    def __init__(self):
        pass
    # INIT

    def init(self,
             loss,
             pdb=None,
             estimator="RF",  # "dummy" for random search
             identifier=None,  # string to identify optimization run
             optimizer=None,    # the optimizer
             distributor=None,  # the distributor
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
             hpc=False,
             signal=None
             ):
        self.identify = identifier
        self.base_estimator = estimator
        self.pandas = save_pandas
        self.loss = loss
        self.pdb = pdb
        self.results = None
        self.signal = signal
        # CONSTANTS
        self.n_cores = n_cores
        self.rpc = rpc
        self.evals = evals
        self.test_run = test_run
        # COUNTER
        self.batch_counter = int(self.n_cores/self.rpc)
        self.evals_done = 0
        # MEMBER CLASSES
        self.distributor = distributor
        self.optimizer = optimizer
        self.logger = logging.getLogger('OptimizationManager')
        self.logger.setLevel(logging.DEBUG)
        # TEST CASE

        if test_run:
            self.logger.info('TEST RUN')
            self.objective = config._dummy_objective.remote
            self.loss_value = 'ref15'
            self.init_method = None
            self.pandas = False
            estimator = 'dummy'
        else:
            self.objective = config._objective
            self.init_method = config._init_method

        self.logger.debug('CHECK SERIALIZABLE ACTOR INITIALIZER FUNCTION')
        inspect_serializability(self.init_method, 'initializer method')
        # SEARCH SPACE

        # OPTIMIZER
        self.optimizer.init.remote(
            random_state=5,
            dimensions=config.space_dimensions,
            base_estimator=estimator,
            acq_func_kwargs=config.acq_func_kwargs,
            n_initial_points=n_cores/rpc,
            cooldown=cooldown,
            evals=evals
        )

        # DISTRIBUTOR
        self.distributor.init.remote(
            manager_callback=self.log_res_and_update,
            hpc=hpc,
            workers=n_cores,
            rpc=rpc,
            evals=evals,
            initializer=self.init_method)

        # BOOKKEEPING
        # make empty array
        # self.results = np.array(self.evals)
        self.logger.debug('initialized Manager')

    def log_res_and_update(self, map_res: list = None, make_batch: bool = True) -> None:
        self.evals_done += 1 # we need separate counters for made and returned batches
        self.results = np.append(self.results, map_res)
        self.logger.debug('self.results %s', self.results)
        if not map_res:
            self.make_batch()
        else:
            if len(map_res) == self.rpc:
                self.optimizer.update_prior.remote(
                    map_res[0]['config'], sum(res[self.loss] for res in map_res)/len(map_res))
            else:
                self.logger.error(
                    'Distributor returned %d results but self.rpc is %d', len(map_res), self.rpc)
                raise Exception('Distributor returned %d results but self.rpc is %d'.format(
                    len(map_res), self.rpc))

            if self.batch_counter < self.evals:
                self.make_batch()
            elif self.evals_done == self.evals:
                self._save_and_exit()

    def make_batch(self):
        config = self.optimizer.get_next_config.remote()
        self.distributor.distribute.remote(func=self.objective,
                                           params=config,
                                           pdb=self.pdb,
                                           num_workers=self.rpc,
                                           run=self.batch_counter+1)
        self.batch_counter += 1

    def _save_and_exit(self) -> bool:
        """
            Saves the results stored in the DataFrame and reports to console
        """
        self.logger.debug('DONE')

        if not self.test_run:
            if self.pandas:
                # save pandas DataFrame with correct column names
                df = pd.DataFrame(self.results.tolist()[1:])
                breakpoint()
                self.logger.debug(df)
                weights = df.config.apply(lambda x: pd.Series(x))
                weights.columns = ["fa_rep_" + str(num) for num in range(7)]
                results = pd.concat([df, weights], axis=1)
            with open(
                config.result_path+"/{}_{}_res_{}.pkl".format(self.identify,
                                                  self.base_estimator, self.evals),
                "wb",
            ) as file:
                pickle.dump(results, file)
        else:
            self.logger.info('len results %d', len(self.results))
        self.logger.warning(
            ';;;;;FINAL STATE;;;;;; \n evals: %s \n  results: %s ', self.evals_done, self.results)
        self.distributor.terminate.remote()
        ray.get(self.signal.send.remote())
        print("TERMINATING")

    def run(self) -> None:
        """
            Runs the Manager and 
        """
        self.logger.debug('RUN OptimizationManager')
        # map initial runs workers/rpc rpc times
        self.logger.debug('CHECK PICKABLE OBJECTIVE')
        inspect_serializability(self.objective, name='objective_func')

        for run in range(int(self.n_cores/self.rpc)):

            self.distributor.distribute.remote(func=self.objective,
                                               params=self.optimizer.get_next_config.remote(),
                                               pdb=self.pdb,
                                               num_workers=self.rpc,
                                               run=run+1)
        self.logger.info('INITIAL DISTRIBUTION DONE GOING TO CYCLIC')
