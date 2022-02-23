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
             pdb=None,     # DEPRECATED: move to args
             target=None,  # DEPRECATED: move to args
             estimator="RF",  # "dummy" for random search
             identifier=None,  # string to identify optimization run
             optimizer=None,    # the optimizer
             distributor=None,  # the distributor
             test_run=False,  # if test run eval dummy_objective instead of real
             evals=200,  # configuration evaluations on the objective
             rpc=8,  # runs_per_config n_calls/rpc = evals
             warm_start=None,  # continue previous optimization run
             n_cores=None,
             cooldown=True,  # cooldown exploration to exploitation
             space_dimensions=None,  # yaml file with optimizer dimensions
             save_pandas=True,  # save results as pickled pandas DataFrames
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
        self.cooldown = cooldown
        # COUNTER
        self.batch_counter = 0
        self.start_time = time.time()
        self.current_cycle = 1
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
            self.base_estimator = 'dummy'
        else:
            self.objective = config._objective
            self.init_method = config._init_method

        self.logger.debug('CHECK SERIALIZABLE ACTOR INITIALIZER FUNCTION')
        inspect_serializability(self.init_method, 'Initializer Method')
        inspect_serializability(self.objective, 'Objective Method')

        # BOOKKEEPING
        self.batches = {}
        self.results = []

        # INITIALIZE
        self.init_distributor()
        self.init_optimizer()

    def init_optimizer(self):
        # OPTIMIZER
        self.optimizer.init(
            random_state=5,
            dimensions=config.space_dimensions,
            base_estimator=self.base_estimator,
            acq_func_kwargs=config.acq_func_kwargs,
            n_initial_points=int(self.n_cores//self.rpc),
            cooldown=self.cooldown,
            evals=self.evals
        )

    def init_distributor(self):
        # DISTRIBUTOR
        self.distributor.init(
            workers=self.n_cores,
            rpc=self.rpc,
            evals=self.evals,
            initializer=self.init_method)

    def no_optimize(self, identify, evals, pdb, config_path=None):
        """RUN objective with (default)config evals times without ommptimization 


        """
        # distribute work evenly over workers
        res = self.distributor.distribute(func=self.objective,
                                          params=None,
                                          pdb=self.pdb,
                                          target=self.target,
                                          run=1,
                                          num_workers=evals,
                                          round_robin=True)
        self.results = self.distributor.get_batch(
            complete_run_batch=False, batch_size=evals)
        self.save()

    def log_res_and_update_optimizer(self, map_res: list = None, make_batch: bool = True) -> None:
        """
        Append the result obtained from the distributor to the results list.
        Calls optimizer to incoorporate the new Info form results.
        """
        self.logger.info('RUN %d DONE', map_res[0]['run'])

        for r in map_res:
            r['cycle'] = self.current_cycle
        self.results.extend(map_res)
        self.optimizer.handle_result(
            map_res[0]['config'], sum(res[self.loss] for res in map_res)/len(map_res))

    def make_batch(self, round_robin=False):
        """makes a new job batch by retrieving a new point x_i form the optimizer
        and calling the distributors distribute function 

        @param round_robin: indicate whether to distribute work round_robin 
        over workers or if False look for idle workers
        @type round_robin: boolean

        @return: None
        @rtype : None

        """
        self.batch_counter += 1
        config = self.optimizer.get_next_config()
        self.distributor.distribute(func=self.objective,
                                    params=config,
                                    pdb=self.pdb,
                                    target=self.target,
                                    num_workers=self.rpc,
                                    run=self.batch_counter,
                                    round_robin=round_robin)

    def get_results(self):
        """
            getter for self.results

        @return: self.results the current state of the result list
        @rtype: pd.DataFrame 
        """

        if self.pandas:
            # print(self.results)
            df = pd.DataFrame(self.results)
            # df['cycle'] = self.current_cycle
            print(df)

            if df.empty:
                raise Exception("Result DataFrame is empty")

            return df
        else:
            return self.results

    def save(self) -> bool:
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
                "{}{}.pkl".format(config.result_path, self.identify),
                "wb",
            ) as file:
                pickle.dump(df, file)
        else:
            self.logger.info('len results %d', len(self.results))
        # TIME DELTA
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

    def add_res_to_batch(self, result, batch_number):
        """adds a single result to a bat

        @param param:  Description
        @type  param:  Type

        @return:  Description
        @rtype :  Type

        @raise e:  Description
        """

        if batch_number in self.batches.keys():
            self.batches[batch_number].append(result)
        else:
            self.batches.update({batch_number: [result]})

    def run(self, report=False, complete_run_batch=False) -> None:
        """
            Runs the Manager and 
        """
        self.logger.info('RUN')
        # map initial runs workers/rpc times

        initial_runs = int(self.n_cores/self.rpc)

        for run in range(initial_runs):
            self.make_batch(round_robin=True)
        self.logger.info('INITIAL DISTRIBUTION DONE GOING TO CYCLIC')

        # Actual RUN LOOP

        for run in range(self.evals):
            # NOTE: blocks until batch is ready and update
            batch = self.distributor.get_batch(
                complete_run_batch=complete_run_batch)

            for r in batch:  # handle result to batch matching
                self.add_res_to_batch(r, r['run'])

                if len(self.batches[r['run']]) == self.rpc:
                    # handle complete batch
                    self.log_res_and_update_optimizer(self.batches[r['run']])

            if run < (self.evals - initial_runs):
                # make batch
                self.make_batch()

        if report:
            return self.get_results()

    def set_pdb(self, pdb):
        self.pdb = pdb

    def set_cycle(self, cycle):
        self.current_cycle = cycle
