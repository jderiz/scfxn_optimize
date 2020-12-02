import logging
import multiprocessing
import os
import pickle
import random
import sys
import time
from functools import partial
from multiprocessing import (Pipe, active_children, cpu_count,
                             current_process, get_context)

from numpy import repeat

from design import design_with_config
from hyperparams import ref15_weights, scfxn_ref15_space
from skopt import (Optimizer, callbacks, forest_minimize, gbrt_minimize,
                   gp_minimize)


class Optimization:
    """
        Optimization class that runs a specified optimizer.
        It Instatiates the Multiprocessing Pool and handles the result callbacks
    """

    # Setup Logging
    logger = multiprocessing.log_to_stderr()
    logger.setLevel(logging.WARNING)

    def __init__(
        self,
        thread_pool,
        loss,
        base_estimator='RF',
        test_run=False,
        n_calls=200,
        runs_per_config=5,
        result_dir="results",
        warm_start=None,
        cores=None,
    ):
        """
        Init Optimization 
        """
        self.tp = thread_pool
        print(self.tp)
        # Setup result folder
        os.makedirs(result_dir, exist_ok=True)
        self.result_buffer = {}
        self.results = []

        # Parallel handling

        if cores:
            self._cores = cores
        else:
            self._cores = cpu_count()
        self.n_calls = n_calls
        self.runs_per_config = runs_per_config
        self.calls = 0  # iteration counter
        self.num_callbacks = 0  # keep track of callback calls
        self.jobs = []
        self._DONE = False
        # overall Runtime measuring
        self.start_time = time.time()
        # Optimizer setup
        self.dimensions = scfxn_ref15_space

        if test_run:
            self.objective = self.dummy_objective
        else:
            self.objective = design_with_config
        self.default_parameters = [val for k, val in ref15_weights]
        self.estimator = base_estimator
        self.loss_value = loss
        # more exploit then explore: DEFAULTS: xi:0.01, kappa:1.96
        # TODO: implement cooldown [start_values, end_values]
        xi = 0.001
        kappa = 0.1
        acq_func_kwargs = {"xi": xi, "kappa": kappa}
        self.optimizer = Optimizer(
            dimensions=self.dimensions,
            base_estimator=base_estimator,
            acq_func_kwargs=acq_func_kwargs,
            # n_initial_points=300
        )
        # Warmstart Optimizer

        if warm_start:
            with open(warm_start + ".pkl", "rb") as h:
                wsres = pickle.load(h)

            y_0 = [x[self.loss_value] for x in wsres]
            x_0 = [x["weights"] for x in wsres]
            self.optimizer.tell(x_0, y_0)

    def dummy_objective(self, config):
        # print('TEST')
        time.sleep(1)

        return {
            "bloss62": random.randint(1, 100),
            "ref15": random.randint(1, 50),
            "scfxn": random.randint(1, 46),
        }

    def save_and_exit(self):
        print("wait for all Processes to finish and close pool")
        print(len(active_children()))

        while not self._DONE:
            for job in self.jobs:
                if not job.ready():
                    time.sleep(60)
                    print(job._job, job._cache, job._event, job._pool)
                    job._event.set()  # tell process in same thread that is is done
                    # print(active_children())

            if all([job.ready() for job in self.jobs]):
                print("ALL JOBS READY")
                self.tp.close()

                break

            if not active_children():
                print("No More active_children")
                self._DONE = True

                break
            time.sleep(5)

        print([job.ready() for job in self.jobs])
        print([(key, len(val)) for key, val in self.result_buffer.items()])
        print(len(self.result_buffer))
        print("SAVING")
        took = time.time() - self.start_time
        print(
            "Took for ALL: {} to run".format(
                time.strftime("%D_%H: %M: %S", time.gmtime(took))
            )
        )
        # save custom results
        with open(
            "results/{}_{}_res_{}_{}.pkl".format(
                time.strftime("%D_%H_%M"),
                self.estimator,
                str(self.n_calls),
                self.loss_value,
            ),
            "wb",
        ) as file:
            pickle.dump(self.results, file)
        print("TERMINATING")
        self.tp.terminate()
        _DONE = True

    def make_process(self):

        # while not reached n_calls make either new config and first job for it
        # of instantiate job for existing config until @param runs_per_config.

        if self.cached_config and self.jobs_for_current_config < self.runs_per_config:
            config = self.cached_config
            self.jobs_for_current_config += 1
        else:  # make new config reset counter to 1
            config = self.optimizer.ask()
            self.optimizer.update_next()
            self.cached_config = config
            self.jobs_for_current_config = 1
        # instantiate job knowing which config it belongs to.
        job = self.tp.apply_async(
            self.objective,
            (config,),
            callback=partial(self._callback, config=config),
            error_callback=self.error_callback,
        )
        # Append map_result to jobs list
        self.jobs.append(job)
        self.calls += 1

    def _callback(self, map_res, config=None):
        """
        Whenever a process finishes it calls config_callback with its result
        and a new process can be run
        """

        self.num_callbacks += 1
        print("CALLBACK: ", self.num_callbacks)
        # only frozenset and tuple are hashable
        # c_hash = hash(frozenset(config))
        # sorted(config)
        c_hash = str(sorted(config))

        # Check if key exists

        if c_hash in self.result_buffer.keys():

            if len(self.result_buffer[c_hash]) == (self.runs_per_config - 1):
                # if this is the last run for config
                print("SAVE RESULT")
                # if so, compute mean over results for config and tell optimizer
                self.result_buffer[c_hash].append(map_res)
                # print(result_buffer[c_hash])
                res = self.result_buffer[c_hash]
                _res = {
                    "bloss62": sum([x["bloss62"] for x in res]) / len(res),
                    "ref15": sum([x["ref15"] for x in res]) / len(res),
                    "scfxn": sum([x["scfxn"] for x in res]) / len(res),
                    "weights": config,
                }
                self.results.append(_res)
                self.optimizer.tell(config, _res[self.loss_value])
                # del result_buffer[c_hash]  # dont need that buffer
                # optimizer_results.append(opti_res)
            else:  # not last run yet, append to buffer for this config
                print("APPEND map_res ")
                self.result_buffer[c_hash].append(map_res)
        else:  # otherwise make new buffer entry
            print("NEW BUFFER")
            self.result_buffer.update({c_hash: [map_res]})

        # while there are calls left on the budget

        if self.calls < self.n_calls:
            # Make New Process
            self.make_process()
        else:
            # wait for all processors to finish then exit Pool
            while active_children():
                print(len(self.jobs))
                time.sleep(60)

            # self.pending_jobs -= 1
            # print("PENDING: ", self.pending_jobs)
            # print("ACTIVE CHILDREN: ", len(active_children()))
            # # print([job.ready() for job in jobs])

            # # TODO: think HARD why there is one left job that is not ready ??
            # # has not returned? counter error?

            # if self.pending_jobs == 0:
            #     self.save_and_exit()

    def error_callback(self, r):
        self.logger.log(level=logging.ERROR, msg=r)

    def start(self):
        """
        start the optimization process, and wait for it to finish
        """
        # maxtasksperchild for controlling memory usage
        self.cached_config = None
        self.jobs_for_current_config = 0

        # initial runs

        for _ in range(self._cores):
            self.make_process()

        while not self._DONE:
            time.sleep(60)
            pass
