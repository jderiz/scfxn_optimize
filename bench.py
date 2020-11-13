import logging
import multiprocessing
import os
import pickle
import sys
import time
from functools import partial
from multiprocessing import (Pipe, Pool, active_children, cpu_count,
                             current_process, get_context)

from joblib import Parallel, delayed
from numpy import repeat

from design import design_with_config
from hyperparams import ref15_weights, scfxn_ref15_space
from setup import parallel_configs, runs_per_config, n_calls
from skopt import (Optimizer, callbacks, forest_minimize, gbrt_minimize,
                   gp_minimize)

logger = multiprocessing.log_to_stderr()
logger.setLevel(logging.WARNING)


if __name__ == "__main__":
    """
    Main execution for optimization run.
    """
    # Setup result folder
    os.makedirs("results", exist_ok=True)

    # instantiate result array and specific number calls to objective per optimizer
    results = []
    # n_calls = 800  # Objective Function evaluations
    calls = 0  # iteration counter
    num_callbacks = 0 # keep track of callback calls
    start_time = time.time()  # overall Runtime measuring
    dimensions = scfxn_ref15_space
    objective = design_with_config
    default_parameters = [val for k, val in ref15_weights]
    # "GP, GBRT, ET, RF"
    estimator = sys.argv[1]
    loss_value = sys.argv[2]

    if not loss_value:
        print('ERROR: No loss value defined as second argument')
        exit() 
    # more exploit then explore
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
    optimizer_results = []
    _DONE = False
    # TODO: spawn or fork processes 
    with get_context('spawn').Pool(processes=cpu_count(), maxtasksperchild=2) as tp:
        # maxtasksperchild for controlling memory usage
        result_buffer = {}
        cached_config = None
        jobs_for_current_config = 0
        # For Checking outstanding jobs befor terminate.
        pending_jobs = cpu_count()
        jobs = []

        def save_and_exit():
            print("wait for all Processes to finish and close pool")
            print(len(active_children()))
            # TODO: be certain all jobs are finished
            # while active_children():
            #     print(len(active_children()))
            #     time.sleep(10)
                # pass
            print('No More active_children')
            print('SAVING')
            took = time.time() - start_time
            print("Took for ALL: {} to run".format(time.strftime("%H: %M: %S", time.gmtime(took))))
            # save custom results
            with open(
                "results/{}_{}_res_{}_{}.pkl".format(
                   time.strftime("%H_%M"), estimator, str(n_calls), loss_value
                   ),
                "wb",
                 ) as file:
                pickle.dump(results, file)
            # save results from optimizer.tell()
            with open(
                "results/{}_{}_res_{}_{}_optimizer.pkl".format(
                   time.strftime("%H_%M"), estimator, str(n_calls), loss_value
                   ),
                "wb",
                ) as file:
                pickle.dump(optimizer_results, file)
            print('TERMINATING')
            tp.terminate()
            tp.close()
            _DONE = True

        def make_process():
            global calls
            global cached_config
            global jobs_for_current_config
            global pending_jobs

            if calls < n_calls:
                # while not reached n_calls make either new config and first job for it
                # of instantiate job for existing config until @param runs_per_config.

                if cached_config and jobs_for_current_config < runs_per_config:
                    config = cached_config
                    jobs_for_current_config += 1
                else:  # make new config reset counter to 0
                    config = optimizer.ask()
                    cached_config = config
                    jobs_for_current_config = 1
                # instantiate job knowing which config it belongs to.
                job = tp.apply_async(
                        objective,
                        (config,),
                        callback=partial(_callback, config=config),
                        error_callback=error_callback,
                        )
                # Append map_result to jobs list
                jobs.append(job)
                calls += 1
            else:
              # wait for all processes to finish then exit Pool
              pending_jobs -= 1
              print('PENDING: ', pending_jobs)

            if pending_jobs == 0:
                save_and_exit()

        def _callback(map_res, config=None):
            """
            Whenever a process finishes it calls config_callback with its result
            and a new process can be run
            """
            global result_buffer, results, num_callbacks

            num_callbacks += 1
            print('CALLBACK: ', num_callbacks)
            # only frozenset and tuple are hashable
            c_hash = hash(frozenset(config))

            # Check if key exists

            if c_hash in result_buffer.keys():
                # if this is the last run for config 

                if len(result_buffer[c_hash]) == (runs_per_config - 1):
                    print('SAVE RESULT')
                    # if so, compute mean over results for config and tell optimizer
                    result_buffer[c_hash].append(map_res)
                    # print(result_buffer[c_hash])
                    res = result_buffer[c_hash]
                    _res = {
                            "bloss62": sum([x["bloss62"] for x in res]) / len(res),
                            "ref15": sum([x["ref15"] for x in res]) / len(res),
                            "scfxn": sum([x["scfxn"] for x in res]) / len(res),
                            "weights": config,
                            }
                    results.append(_res)
                    opti_res = optimizer.tell(config, _res[loss_value])
                    # del result_buffer[c_hash]  # dont need that buffer
                    optimizer_results.append(opti_res)
                else: # not last run yet, append to buffer for this config
                    print('APPEND map_res ')
                    result_buffer[c_hash].append(map_res)
            else:
                print('NEW BUFFER')
                result_buffer.update({c_hash: [map_res]})
            # Make New Process
            make_process()


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

        # Start Loop

        for _ in range(cpu_count()):
        # map runs_per_config processes with same config
            make_process()

        # while there are jobs in the list being processed wait for them to finish
        # TODO: Understand this better, and maybe refactor

        while not _DONE:
            if not active_children():
                print('NO MORE ACTIVE CHILDREN')
                # when no more child processes sleep a minute to let things settle
                # then terminate Pool
                print([(key, len(val)) for key, val in result_buffer.items()])
                time.sleep(15)
                tp.terminate()
                break
            # for async_result in jobs:
            #     if async_result.ready():
            #         jobs.remove(async_result)

    ##########################
    # result printing and saving
    # took = time.time() - start_time
    # print("Took for ALL: {} to run".format(time.strftime("%H: %M: %S", time.gmtime(took))))
    # # save custom results
    # with open(
    #     "results/{}_{}_res_{}_{}.pkl".format(
    #        time.strftime("%H_%M"), estimator, str(n_calls), loss_value
    #        ),
    #     "wb",
    #      ) as file:
    #     pickle.dump(results, file)
    # # save results from optimizer.tell()
    # with open(
    #     "results/{}_{}_res_{}_{}_optimizer.pkl".format(
    #        time.strftime("%H_%M"), estimator, str(n_calls), loss_value
    #        ),
    #     "wb",
    #     ) as file:
    #     pickle.dump(optimizer_results, file)
