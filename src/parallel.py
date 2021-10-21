import logging
from functools import partial

import ray
from ray.util.multiprocessing import Pool

# from multiprocessing import Pool, get_context



@ray.remote
class Distributor():

    """
    This class manages the parallelization.
    """

    def __init__(self):
        pass

    def init(self, manager_callback, cpus, hpc, initializer):
        self.hpc = hpc

        self.result_list = []
        self.batches = {}
        self.batch_number = 0
        self.manager_callback = manager_callback

        # self.mp = get_context('spawn').Pool(
        # workers, initializer=initializer)
        # self.mp = Pool(processes=cpus)

        self.logger = logging.getLogger('Distributor')
        self.logger.setLevel(logging.DEBUG)
    def distribute(self, func, params, pdb, num_workers, run):
        """
            distribute a function 
        """
        # jobs = self.mp.starmap_async(
        #     func,
        #     [(params, run, pdb)for _ in [0]*num_workers],
        #     callback=self._callback,
        #     error_callback=self._error_callback,
        # )
        batch_refs = [func.remote(params, run, pdb) for _ in range(num_workers)]
        # Append AsyncResult to jobs list
        # self.logger.debug('%s', jobs)
        #TODO: make into list dann pop und put oder sowas
        self.batches.update({run: batch_refs})
        self.batch_number += 1

    def _callback(self, result):
        # logger.debug("%s", ray.get(result[0]))
        self.manager_callback(map_res=result)
        self.result_list.append(result)

    def _error_callback(self, msg):
        self.logger.error(msg)
        exit(msg)

    def terminate(self):
        self.logger.debug('TERMINATE Pool')
        # self.mp.close()
        # self.mp.terminate()

    def has_batch(self) -> bool:
        """
        returns true if there is something in the result list
        """

        return not self.result_list

    def get_result(self) -> dict:
        '''
        pop result batch from queue and return to manager
        '''

        return self.result_list.pop(0)
    
    def get_result_list(self) -> list:
        return self.result_list

    def get_jobs(self) -> list:
        '''
        return the list of jobs
        '''

        return self.jobs
