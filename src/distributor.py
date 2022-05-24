import logging
import os
import random

import ray
from ray.util import inspect_serializability

from rosetta_worker import PRSActor


class Distributor():

    """
    This Class manages the parallelization. 
    """

    def __init__(self):
        pass

    def init(self, workers, rpc, evals, init_args=None):
        self.logger = logging.getLogger('Distributor')
        self.logger.setLevel(logging.DEBUG)
        # init batch dict with empty list for each batch
        self.batches = {}  # storing single results
        self.batch_size = rpc
        self.rpc_counter = 0
        self.futures = []  # this holds our object_refs to the worker tasks
        self.run_futures = {}
        # MAKE POOL
        self._initargs = init_args
        self._start_actor_pool(workers)
        self.next_rr_idx = 0
        # self.mp = Pool(processes=workers, initializer=initializer)
        self.logger.debug('CWD %s', os.getcwd())
        self.logger.info('INTITIALIZED DISTRIBUTOR')

    def _start_actor_pool(self, processes):
        # make all but one actor have separate ressources and one actor that shares with framework
        self._actor_pool = [self._new_actor_entry(
            num_cpus=1, idx=idx) for idx in range(processes-1)]
        self._actor_pool.append(self._new_actor_entry(False, idx=processes-1))
        ray.get([actor.ping.remote() for actor, _ in self._actor_pool])

    def _new_actor_entry(self, num_cpus, idx=None):
        # NOTE(edoakes): The initializer function can't currently be used to
        # modify the global namespace (e.g., import packages or set globals)
        # due to a limitation in cloudpickle.
        # HACK: all actors but one that shares with manager, distributor etc
        # are exclusive on cpu

        if num_cpus:
            return(PRSActor.options(num_cpus=1).remote(self._initargs, idx=idx), 0)
        else:
            return (PRSActor.remote(self._initargs, idx=idx), 0)

    def _random_actor_index(self):
        return random.randrange(len(self._actor_pool))

    def _round_robin_index_next(self):
        self.next_rr_idx += 1

        return self.next_rr_idx % len(self._actor_pool)

    def _idle_actor_index(self):
        # self.logger.debug('Get IDLE Actor')
        try:
            # found idle actor, return its index
            idx_ready, _ = ray.wait([actor.ping.remote() for actor, _ in self._actor_pool],
                                    num_returns=1, timeout=5)

            if not idx_ready:
                return None

            return ray.get(idx_ready[0])
        except ray.exceptions.GetTimeoutError:
            return None  # found no idle actor

    def evaluate_config(self, config, fargs, run, error_callback=None,
                        callback=None, round_robin=False) -> tuple:

        if round_robin:
            actor_idx = self._round_robin_index_next()
        else:
            actor_idx = self._idle_actor_index()

            if not actor_idx:  # get random Actor if no idle found
                actor_idx = self._round_robin_index_next()

        if actor_idx != None:
            self.logger.debug("Eval on Actor %d", actor_idx)
            actor, count = self._actor_pool[actor_idx]
            object_ref = actor.evaluate_config.remote(
                config, fargs, run)

            return object_ref
        else:
            raise Exception('could not find an actor_idx to use')

    def distribute(self, config, fargs, run, num_workers=None, round_robin=False):
        """
            distribute a function to the Pool and hold the object_refs 
            to the results somewhere such that done, 
            undone = ray.wait(all_refs, num_returns=batch_size) --> ray.get(the_done_ones) retrieves the results 
        """

        if not num_workers:
            num_workers = self.batch_size
        run_futures = []

        for _ in [0]*num_workers:
            # distrbute to workers
            object_ref = self.evaluate_config(
                config=config,
                fargs=fargs,
                run=run,
                error_callback=self._error_callback,
                round_robin=round_robin,
            )
            # save object_ref to the future list
            self.futures.append(object_ref)
            run_futures.append(object_ref)
        # save the results to the batch dict
        self.run_futures.update({run: run_futures})

    def add_res_to_batch(self, result, batch_number):
        """
        add a result to the batch dict
        @param result: the result to add
        @param batch_number: the batch number to add to
        """

        if batch_number in self.batches.keys():
            self.batches[batch_number].append(result)
        else:
            self.batches.update({batch_number: [result]})

    def get_batch(self, complete_run_batch, batch_size=None) -> list:
        """
        block until batch_size tasks are completed and return their result. 
        the rest of the futures becomes the new futures list where new tasks 
        get appended to.

        @param complete_run_batch: if True then the first completed configuration 
        is returned as batch
        @type round_robin: boolean

        @param batch_size: if not complete_run_batch indicate how many results
        should be included in the batch. usually runs_per_config many
        @type batch_size: int

        @return: result batch
        @rtype: list

        """
        # if no batch size is given, use the default batch size
        batch_size = batch_size if batch_size else self.batch_size

        if complete_run_batch:
            # blocks until an entire run_batch is complete
            # only sensible when runs take aprox. the same time.
            self.logger.debug('Collect complete_run_batch')
            wait = True

            while wait:  # wait a run is complete, NOTE: this is a blocking call
                for run_batch_key in self.run_futures.keys():
                    # try the run futures one after another
                    ready_batch, undone = ray.wait(
                        self.run_futures[run_batch_key],
                        timeout=3, num_returns=batch_size)

                    if len(ready_batch) == batch_size:  # found complete batch
                        del self.run_futures[run_batch_key]  # delete from dict
                        wait = False  # exit while loop, as complete run batch found

                        break  # exit for-loop
                    else:
                        continue
        else:
            # blocks until batch_size results are done
            ready_batch, undone = ray.wait(
                self.futures, num_returns=batch_size)
            self.futures = undone  # set futures to the remaining
        self.logger.debug('READY: %s', ready_batch)
        self.logger.debug('UNDONE: %s', undone)
        batch = ray.get(ready_batch)
        self.logger.debug(batch)

        return batch

    def _error_callback(self, msg):
        self.logger.error(msg)

    def terminate_pool(self):
        for a, c in self._actor_pool:
            a.shutdown.remote()
