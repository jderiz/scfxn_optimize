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
from bayesopt import BayesOpt
from distributor import Distributor


# @ray.remote
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
             target=None,
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
             cooldown=True,  # cooldown exploration to exploitation
             space_dimensions=None,  # yaml file with optimizer dimensions
             save_pandas=True, # save results as pickled pandas DataFrames
             ):
        self.identify = identifier
        self.base_estimator = estimator
        self.pandas = save_pandas
        self.loss = loss
        self.pdb = pdb
        self.target = target
        self.results = None
        # CONSTANTS
        self.n_cores = n_cores
        self.rpc = rpc
        self.evals = evals
        self.test_run = test_run
        # COUNTER
        self.batch_counter = 0
        self.start_time = time.time()
        # MEMBER CLASSES
        self.distributor: Distributor = distributor
        self.optimizer: BayesOpt = optimizer

        self.logger = logging.getLogger('OptimizationManager')
        self.logger.setLevel(logging.DEBUG)
        self.logger.debug('CWD %s', os.getcwd())
        # TEST CASE

        if test_run:
            self.logger.info('TEST RUN')
            self.objective = config._dummy_objective
            config._init_method = None
            config._objective = config._dummy_objective
            self.loss_value = 'score'
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
        self.optimizer.init(
            random_state=5,
            dimensions=config.space_dimensions,
            base_estimator=estimator,
            acq_func_kwargs=config.acq_func_kwargs,
            n_initial_points=int(n_cores//rpc),
            cooldown=cooldown,
            evals=evals
        )

        # DISTRIBUTOR
        self.distributor.init(
            manager_callback=self.log_res_and_update,
            workers=n_cores,
            rpc=rpc,
            evals=evals,
            initializer=self.init_method)

        # BOOKKEEPING
        self.batches = {}
        self.results = []

    def no_optimize(self, identify, evals, pdb, config_path=None):
        """RUN objective with (default)config evals times without ommptimization 

        :arg1: TODO
        :returns: TODO

        """
        res = self.distributor.distribute(func=self.objective, params=None,
                                          pdb=pdb, run=1, num_workers=evals, round_robin=True)
        self.results = self.distributor.get_batch(batch_size=evals)
        self._save_and_exit()

    def log_res_and_update(self, map_res: list = None, make_batch: bool = True) -> None:
        self.logger.info('RUN %d DONE', map_res[0]['run'])
        self.results.extend(map_res)
        self.optimizer.handle_result(
            map_res[0]['config'], sum(res[self.loss] for res in map_res)/len(map_res))

    def make_batch(self, round_robin=False):
        self.batch_counter += 1
        config = self.optimizer.get_next_config()
        self.distributor.distribute(func=self.objective,
                                    params=config,
                                    pdb=self.pdb,
                                    target=self.target,
                                    num_workers=self.rpc,
                                    run=self.batch_counter,
                                    round_robin=round_robin)

    def report(self):
        """
        reports the results
        :returns: self.results
        """

        if self.pandas:
            return pd.DataFrame(self.results)
        else:
            return self.results

    def _save_and_exit(self) -> bool:
        """
            Saves the results stored in the DataFrame and reports to console
        """
        self.logger.debug('DONE \n test %s \n results %d',
                          self.test_run, len(self.results))

        if not self.test_run:
            if self.pandas:
                # save pandas DataFrame with correct column names
                df = pd.DataFrame(self.results)
                self.logger.debug(df)
            with open(
                config.result_path+"/{}_{}_res_{}.pkl".format(self.identify,
                                                              self.base_estimator, self.evals),
                "wb",
            ) as file:
                pickle.dump(df, file)
        else:
            self.logger.info('len results %d', len(self.results))
        tdelt = time.time()-self.start_time
        days = tdelt//86400
        tdelt = tdelt - (days*86400)
        hours = tdelt//3600
        tdelt = tdelt - (hours*3600)
        minutes = tdelt//60
        took = "{} days {} hours {} minutes".format(days, hours, minutes)
        self.logger.info(
            '\n -------- FINAL STATE -------- \n Got %d Results \n TOOK %s',
            len(self.results), took)
        self.distributor.report()
        self.optimizer.report()
        self.distributor.terminate_pool()

    def add_res_to_batch(self, result, batch_number):
        if batch_number in self.batches.keys():
            self.batches[batch_number].append(result)
        else:
            self.batches.update({batch_number: [result]})

    def run(self, report=False) -> None:
        """
            Runs the Manager and 
        """
        self.logger.info('RUN')
        # map initial runs workers/rpc rpc times

        initial_runs = int(self.n_cores/self.rpc)

        for run in range(initial_runs):
            self.make_batch(round_robin=True)
        self.logger.info('INITIAL DISTRIBUTION DONE GOING TO CYCLIC')
        # REfact

        for run in range(self.evals):
            # blocks until batch is ready and update
            batch = self.distributor.get_batch()

            for r in batch:
                self.add_res_to_batch(r, r['run'])

                if len(self.batches[r['run']]) == self.rpc:
                    self.log_res_and_update(self.batches[r['run']])

            # make batch

            if run < (self.evals - initial_runs):
                self.make_batch()
        # finally report optimization result or save the result and exit

        if report:
            return self.report()
        self._save_and_exit()
