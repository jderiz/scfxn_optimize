import logging
from functools import partial

import ray
from ray.util import inspect_serializability

from mypool import Pool

# from rosetta_worker import PRSActor


# @ray.remote
class Distributor():

    """
    This Actor manages the parallelization.
    """

    def __init__(self):
        pass

    def init(self, manager_callback, workers, rpc, evals, hpc, initializer):
        self.logger = logging.getLogger('Distributor')
        self.logger.setLevel(logging.DEBUG)
        # init batch dict with empty list for each batch
        self.batches = {}  # storing single results
        self.manager_callback = manager_callback
        self.batch_size = rpc
        self.logger.info('INTITIALIZED DISTRIBUTOR')
        self.logger.info(' batches %s', self.batches)
        self.rpc_counter = 0
        self.futures = []  # this holds our object_refs to the worker tasks
        # MAKE POOL
        self.mp = Pool(processes=workers, initializer=initializer)

    def distribute(self, func, params, pdb, run, num_workers=None):
        """
            distribute a function to the Pool and hold the object_refs to the results somewhere such that done, undone = ray.wait(all_refs, num_returns=batch_size) --> ray.get(the_done_ones) retrieves the results 
        """

        if not num_workers:
            num_workers = self.batch_size
        # self.logger.debug('map config for run %s %s times', run, num_workers)

        for _ in [0]*num_workers:
            # self.logger.debug("run %d apply_async %d", run, _)
            object_ref = self.mp.evaluate_config(
                params=params,
                run=run,
                pdb=pdb,
                error_callback=self._error_callback,
            )
            self.futures.append(object_ref)
        self.logger.debug(len(self.futures))

    def add_res_to_batch(self, result, batch_number):
        if batch_number in self.batches.keys():
            self.batches[batch_number].append(result)
        else:
            self.batches.update({batch_number: [result]})

    def get_batch(self):
        """
        block until self.batches tasks are completed and return their result. the rest of the futures becomes the new futures list where new tasks get appended to.
        """
        # blocks until batch_size results are done 
        ready_batch, undone = ray.wait(
            self.futures, num_returns=self.batch_size)
        self.logger.debug(ready_batch)
        batch = ray.get(ready_batch)
        self.logger.debug(batch)
        self.futures = undone

        return batch

    def _callback(self, result):
        """
        Aggregate batch results and call manager callback once rpc workers have returned a result, 
        when a batch is completed distinguish if we just update the BayesOpt or if we also make a new batch
        """

        self.rpc_counter += 1
        run = result['run']
        self.add_res_to_batch(result, run)
        try:
            self.logger.warning('%s', [(key, len(val))
                                       for key, val in self.batches.items()])
        except Exception as e:
            self.logger.error('len self.batches %s %s',
                              len(self.batches), self.batches)

        if (len(self.batches[run]) == self.batch_size):
            # batch complete

            if (self.rpc_counter == self.batch_size):
                # call manager_callback with completed batch and remove batch from batches
                self.logger.debug('CALLING MANAGER CALLBACK for run %s', run)
                self.manager_callback(map_res=self.batches[run])
                self.logger.debug(
                    'removing key %s from batches dict', run)
                del self.batches[run]
                self.rpc_counter = 0
            else:  # call manager callback to update BayesOpt but dont make a new batch just jet.
                self.manager_callback(
                    map_res=self.batches[run], make_batch=False)
                del self.batches[run]
        elif self.rpc_counter == self.batch_size:
            # make new batch but dont update optimizer
            self.manager_callback(map_res=None)
            self.rpc_counter = 0

    def _error_callback(self, msg):
        self.logger.error(msg)

    def terminate(self):
        self.logger.debug('TERMINATE Pool')
        self.mp.close()
        self.mp.terminate()
