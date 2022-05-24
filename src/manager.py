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

import config
from bayesopt import BayesOpt
from distributor import Distributor


class OptimizationManager():
    """
    Optimization Manager handles the communication between optimizer and 
    distributor sets up the optimization and handles the results.
    """
    # dummy init

    def __init__(self):
        pass
    # INIT

    def init(self,
             loss,
             fargs=None,
             estimator="RF",  # "dummy" for random search
             identifier=None,  # string to identify optimization run
             optimizer=None,    # the optimizer
             distributor=None,  # the distributor
             test_run=False,  # if test run eval dummy_objective instead of real
             evals=200,  # configuration evaluations on the objective
             rpc=6,  # runs_per_config n_calls/rpc = evals
             warm_start=None,  # continue previous optimization run
             n_cores=None,
             cooldown=True,  # cooldown exploration to exploitation
             space_dimensions=None,  # yaml file with optimizer dimensions
             save_pandas=True,  # save results as pickled pandas DataFrames
             ):
        self.fargs = fargs
        self.identify = identifier
        self.base_estimator = estimator
        self.pandas = save_pandas
        self.loss = loss
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

        # BOOKKEEPING
        self.batches = {}
        self.results = []

        # INITIALIZE
        self.init_distributor()
        self.init_optimizer()

    def init_optimizer(self):
        # OPTIMIZER instantiation
        self.optimizer.init(
            random_state=5,
            base_estimator=self.base_estimator,
            acq_func_kwargs=config.acq_func_kwargs,
            n_initial_points=int(self.n_cores // self.rpc),
            cooldown=self.cooldown,
            evals=self.evals
        )

    def init_distributor(self):
        # DISTRIBUTOR instantiation
        self.distributor.init(
            workers=self.n_cores,
            rpc=self.rpc,
            evals=self.evals,
        )

    def no_optimize(self, evals, run, config_path=None):
        """
        RUN objective with (default)config evals times without ommptimization. 
        @param evals: number of evaluations
        @type evals: int
        @param run: run number
        @type run: int
        @param config_path: path to config file
        @type config_path: str
        """

        if config_path is not None:
            config = pickle.load(open(config_path, 'rb'))
        else:
            config = None
        # distribute work evenly over workers
        res = self.distributor.distribute(
            config=config,
            fargs=self.fargs,
            run=run,
            num_workers=evals,
            round_robin=True)
        result = self.distributor.get_batch(
            complete_run_batch=False, batch_size=evals)
        self._add_result(result)

    def log_res_and_update_optimizer(self, map_res: list = None, make_batch: bool = True) -> None:
        """
        Append the result obtained from the distributor to the results list.
        Calls optimizer to incoorporate the new Info form results.
        @param map_res: list of results
        @type map_res: list
        @param make_batch: indicate whether to make a new batch
        @type make_batch: boolean
        """
        self.logger.info('RUN %d DONE', map_res[0]['run'])

        if map_res is not None:
            self._add_result(map_res)
            self.update_optimizer(map_res)
        elif map_res is None:
            raise Exception('No Result was supplied')

    def _add_result(self, map_res: list) -> None:
        """ 
        Add a result to the results list.
        Appeding the current cycle to the result.
        @param map_res: list of results
        @type map_res: list
        """

        for r in map_res:
            r['cycle'] = self.current_cycle
        self.results.extend(map_res)

    def update_optimizer(self, map_res):
        """
        Update the optimizer with the results of the current batch.
        @param map_res: list of results
        @type map_res: list
        """
        self.optimizer.handle_result(
            map_res[0]['config'], sum(res[self.loss] for res in map_res)/len(map_res))

    def make_batch(self, round_robin=False):
        """makes a new job batch by retrieving a new configuration c_i form the optimizer
        and calling the distributors distribute function

        @param round_robin: indicate whether to distribute work round_robin 
        over workers or if False look for idle workers
        @type round_robin: boolean

        @return: None
        @rtype : None

        """
        self.batch_counter += 1
        config = self.optimizer.get_next_config()
        self.distributor.distribute(
            config=config,
            fargs=self.fargs,
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
            df = pd.DataFrame(self.results)
            # df['cycle'] = self.current_cycle
            self.logger.debug(df)

            if df.empty:
                raise Exception("Result DataFrame is empty")

            return df
        else:
            return self.results

    def save(self) -> bool:
        """
            Saves the results stored in the DataFrame and reports to console
        """
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

        self.logger.debug('DONE \n test %s \n results %d',
                          self.test_run, len(self.results))
        result = self.get_results()
        with open(
            "{}{}.pkl".format(config.result_path, self.identify),
            "wb",
        ) as file:  # save all except pose
            pickle.dump(result.drop('pose', axis=1), file)

        with open(
            "{}{}.pkl".format(config.pose_dir, self.identify),
            "wb",
        ) as file:  # save pose
            pickle.dump(result.pose, file)

    def add_res_to_batch(self, result, batch_number):
        """adds a single result to a batch.
        @param result: the result to be added
        @type result: dict
        @param batch_number: the batch number to which the result should be added
        @type batch_number: int
        """

        if batch_number in self.batches.keys():
            self.batches[batch_number].append(result)
        else:
            self.batches.update({batch_number: [result]})

    def run(self, report=False, complete_run_batch=False) -> None:
        """
            Runs the Optimization procedure. This function blocks until
            optimization is complete.
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

    def set_fargs(self, fargs):
        self.fargs = fargs

    def set_cycle(self, cycle):
        self.current_cycle = cycle
