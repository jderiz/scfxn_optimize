import logging
from functools import partial
# from multiprocessing import Pool, get_context

from ray.util.multiprocessing import Pool
import ray

# print(ray.nodes(), ray.cluster_resources())

logger = logging.getLogger('Distributor')
logger.setLevel(logging.DEBUG)

# @ray.remote
class Distributor():

    """
    This class manages the parallelization.
    """

    def __init__(self, manager_callback, hpc, cpus, initializer):
        self.hpc = hpc

        self.result_list = []
        self.jobs = []
        self.batch_number = 0
        self.manager_callback = manager_callback

        # self.mp = get_context('spawn').Pool(
            # workers, initializer=initializer)
        self.mp = Pool(cpus, ray_address='auto')

    def distribute(self, func, params, pdb, num_workers, run):
        job = self.mp.starmap_async(
            func,
            [(params, run, pdb)for _ in range(num_workers)],
            callback=self._callback,
            error_callback=self._error_callback,
        )
        # Append map_result to jobs list
        self.jobs.append(job)
        self.batch_number += 1

    def _callback(self, result):
        self.manager_callback(map_res=result)
        self.result_list.append(result)

    def _error_callback(self, msg):
        logger.error(msg)
        exit(msg)

    def terminate(self):
        logger.debug('TERMINATE Pool')
        self.mp.close()
        self.mp.terminate()

    def has_result(self) -> bool:
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
