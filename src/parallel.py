import logging
from functools import partial

import ray
from mypool import Pool
from ray.util import inspect_serializability
# from rosetta_worker import PRSActor

@ray.remote
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
        # minus our three actors dist, manage, opti
        # self.workers = [PRSActor.remote() for _ in range(workers)]
        self.mp = Pool(processes=workers, initializer=initializer)
        self.logger.info('INTITIALIZED DISTRIBUTOR')
        self.logger.info(' batches %s', self.batches)

    def distribute(self, func, params, pdb, run, num_workers=None):
        """
            distribute a function 
        """

        if not num_workers:
            num_workers = self.batch_size
        self.logger.debug('map config for run %s %s times', run, num_workers)
        for _ in [0]*num_workers:
            self.mp.apply_async(
                func,
                args=(params, run, pdb),
                # [(params, run, pdb)for _ in [0]*num_workers],
                callback=self._callback,
                error_callback=self._error_callback,
            )

    def add_res_to_batch(self, result, batch_number):
        if batch_number in self.batches.keys():
            self.batches[batch_number].append(result)
        else:
            self.batches.update({batch_number: [result]})

    def _callback(self, result):
        """
        Aggregate batch results and call manager callback once a batch is complete
        """

        # result = ray.get(result)

        if isinstance(result, list):
            # handle properly this callback is only called once all jobs finish from pool batch
            self.manager_callback(result)

        else:
            # save single result in batch and call manager_callback once #rpc results
            run = result['run']

            self.add_res_to_batch(result, run)
            try:
                self.logger.warning('%s', [(key, len(val))
                                           for key, val in self.batches.items()])
            except Exception as e:
                self.logger.error('len self.batches %s %s',
                                  len(self.batches), self.batches)

            if len(self.batches[run]) == self.batch_size:
                # call manager_callback with completed batch and remove batch from batches
                self.logger.info('CALLING MANAGER CALLBACK for run %s', run)
                self.manager_callback(map_res=self.batches[run])
                self.logger.debug(
                    'removing key %s from batches dict', run)
                del self.batches[run]

    def _error_callback(self, msg):
        self.logger.error(msg)

    def terminate(self):
        self.logger.debug('TERMINATE Pool')
        self.mp.close()
        self.mp.terminate()
