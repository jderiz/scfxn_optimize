import logging
import multiprocessing
import os
import pickle
import sys
import time
from functools import partial
from multiprocessing import Pool, cpu_count, get_context

from joblib import Parallel, delayed
from numpy import repeat
from skopt import (Optimizer, callbacks, forest_minimize, gbrt_minimize,
                   gp_minimize)

from design import design_with_config
from hyperparams import ref15_weights, scfxn_ref15_space
from setup import parallel_configs, runs_per_config

logger = multiprocessing.log_to_stderr()
logger.setLevel(logging.INFO)
logger.warning('doomed')

if __name__ == "__main__":
    """
    Main execution for optimization run.
    """
    # Setup result folder
    os.makedirs("results", exist_ok=True)

    # instantiate result array and specific number calls to objective per optimizer
    res = []
    n_calls = 50  # Objective Function evaluations
    calls = 0  # iteration counter
    start_time = time.time()  # overall Runtime measuring
    dimensions = scfxn_ref15_space
    objective = design_with_config
    default_parameters = [val for k, val in ref15_weights]
    # setup callbacks for logging
    timer_callback = callbacks.TimerCallback()
    forest_check = callbacks.CheckpointSaver(".forest_checkpoints.gz")
    gbrt_check = callbacks.CheckpointSaver(".gbrt_checkpoints.gz")
    gp_check = callbacks.CheckpointSaver(".gp_checkpoints.gz")
    # "GP, GBRT, ET, RF"
    estimator = sys.argv[1]
    xi = 0.001
    kappa = 0.1
    acq_func_kwargs = {"xi": xi, "kappa": kappa}
    print(
        "_________start optimize________"
        + "_____________{}________________".format(estimator)
    )

    optimizer = Optimizer(
        dimensions=dimensions,
        base_estimator=estimator,
        n_jobs=-1,
        acq_func_kwargs=acq_func_kwargs,
    )

    # MAIN optimiaztion Loop with custom result handling
    # using all available cores with Pool()
    with get_context("spawn").Pool(processes=cpu_count()) as tp:
        result_buffer = {}
        cached_config = None
        jobs_for_current_config = 0

        print("initialized Worker Pool with max CPUS")
        jobs = []

        def join_and_exit():
            print("wait for all Processes to finish and close pool")
            tp.join()
            tp.close()

        def make_config_batch():
            global calls
            global cached_config
            global jobs_for_current_config
            print('Making New Batch \n \n \n <>>>>><<<<>')
            print('CALLS TO GO: ', (n_calls - calls))

            # while not reached n_calls

            if calls <= n_calls:
                # check if new config or cached

                if cached_config and jobs_for_current_config <= 3:
                    config = cached_config
                    jobs_for_current_config += 1
                else:  # make new config reset counter to 1
                    config = optimizer.ask()
                    cached_config = config
                    jobs_for_current_config = 1
                # instantiate job with config
                r = tp.apply_async(
                    objective,
                    (config,),
                    callback=partial(config_callback, config=config),
                    error_callback=error_callback,
                )
                # print('Start r.wait()')
                # r.wait(timeout=5.0)
                # results.append(r)
                # Append map_result to jobs list
                jobs.append(r)
                calls += 1
                print('increment calls')

                # return 0
            else:
                join_and_exit()
            print('END of make_config_batch')

        def config_callback(res, config=None):
            """
            Whenever a process finishes it calls config_callback with its result
            and a new process can be run
            """
            print("map_async returning")
            # only frozenset and tuple are hashable
            c_hash = hash(set(config))
            print(c_hash)
            # check if last result for config

            if len(result_buffer[c_hash]) == 3:
                # if so, compute mean over results for config and tell optimizer
                result_buffer[c_hash].append(res)
                _res = {
                    "bloss62": sum([x["bloss62"] for x in res]) / len(res),
                    "ref15": sum([x["ref15"] for x in res]) / len(res),
                    "scfxn": sum([x["scfxn"] for x in res]) / len(res),
                    "weights": config,
                }
                # print(_res)
                res.append(_res)
                print('\n \n <><><<<<<<><> \n CONFIG:', config, '\n RES:', res)
                r = optimizer.tell(config, _res["bloss62"])
                print(r)
            # optimizer.update_next()  # check if necessary
            # make new batch of design processes for next config
            make_config_batch()

        def error_callback(r):
            logger.log(level=logging.ERROR, msg=r)

        # TODO: check for rounding
        initial_batches = int(cpu_count() / runs_per_config)
        print(
            "Starting ",
            initial_batches,
            "initial_configs each running ",
            runs_per_config,
        )

        initial_map_results = []
        prev_config = None
        for _ in range(initial_batches):
            # map runs_per_config processes with same config calling callback\
            # once each is done so a new process can get started
            for _ in range(runs_per_config):
                config = optimizer.ask()

                if config == prev_config:
                    print('SAME CONFIG TWICE:: TERMINATE')
                    exit()
                r = tp.apply_async(
                    objective,
                    (config,),
                    callback=partial(config_callback, config=config),
                    error_callback=error_callback,
                )
                jobs.append(r)
            # time.sleep(1)
            print('\n \n INITIAL :', _)
            print(config)

        # for r in initial_map_results:
        #     r.wait()

        while jobs:
            print('waiting for all jobs to FINISH')

            for r in jobs:
                print(len(jobs))
                r.wait()
                # time.sleep(5)
                jobs.remove(r)
                print(len(jobs))

    # result printing and saving
    took = time.time() - start_time
    print("Took: {} to run".format(time.strftime("%H: %M: %S", time.gmtime(took))))
    with open(
        "results/{}_{}_res_{}.pkl".format(
            time.strftime("%H_%M"), estimator, str(n_calls)
        ),
        "wb",
    ) as file:
        pickle.dump(res, file)
