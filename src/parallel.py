from multiprocessing import (Pipe, Pool, active_children, cpu_count,
                             current_process, get_context)

from dask_jobqueue import SLURMCluster


class Distributor():
    """This class manages the parallelization of work and reports back to the manager"""

    def __init__(self, hpc):
        if not hpc:
            self.setup_single_node()
        else:
            self.setup_multi_node()
        self.result_list = []
        self.jobs = []

    def setup_single_node(n_cores):
        self.mp = get_context('spawn').Pool(
            n_cores, initializer=init_method, maxtasksperchild=mtpc)

    def setup_multi_node(n_cores):
        self.mp = SLURMCluster(
            core=n_cores, mem_per_cpu="6G", partition=partition)
        pass

    def distribute(func, params, num_workers):
        if hpc:

            job = self.mp.map_async(
                partial(func, prot_name),
                [config for _ in range(num_workers)],
                callback=partial(_callback, config=params, run=calls),
                error_callback=error_callback,
            )
            # Append map_result to jobs list
            jobs.append(job)

    def _callback(result):
        self.result_list.append(result)
        pass

    def has_result() -> bool:
        """
        returns true if there is something in the result list
        """

        return not self.result_list

    def get_result() -> dict:
        return self.result_list.pop(0)

    def get_jobs() -> list:
        return self.jobs
