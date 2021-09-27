from multiprocessing import (Pipe, Pool, active_children, cpu_count,
                             current_process, get_context)

from dask_jobqueue import SLURMCluster
from config import _init_method 
from functools import partial
import logging 



class Distributor():

    """
    This class manages the parallelization of work and 
    holds results in a queue until the manager fetches them
    """

    def __init__(self, callback, hpc, workers):
        self.hpc = hpc
        if not hpc:
            self.setup_single_node(workers)
        else:
            self.setup_multi_node(workers)
        self.result_list = []
        self.jobs = []
        self.batch_number = 0
        self.logger = logging.getLogger('Distributor')
        self.manager_callback = callback

    def setup_single_node(self, workers):
        self.mp = get_context('spawn').Pool(
            workers, initializer=_init_method)

    def setup_multi_node(self, workers):
        # TODO: implement
        self.mp = SLURMCluster(
            core=workers, mem_per_cpu="6G")
        pass

    def distribute(self, func, params, num_workers, run):

        job = self.mp.map_async(
            func,
            [params for _ in range(num_workers)],
            callback=partial(self._callback, config=params, run=run),
            error_callback=self._error_callback,
        )
        # Append map_result to jobs list
        self.jobs.append(job)
        self.batch_number +=1

    def _callback(self, config, run, result):
        self.manager_callback(result)
        self.result_list.append(result)
       
    def _error_callback(self, msg):
        self.logger.log(msg)
        print('ERROR: '+ msg)
        exit(msg)


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

    def get_jobs(self) -> list:
        '''
        return the list of jobs 
        '''
        return self.jobs
