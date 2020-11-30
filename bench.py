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
from skopt import (Optimizer, callbacks, forest_minimize, gbrt_minimize,
                   gp_minimize)

from design import design_with_config
from hyperparams import ref15_weights, scfxn_ref15_space
from setup import n_calls, parallel_configs, runs_per_config

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
    num_callbacks = 0  # keep track of callback calls
    start_time = time.time()  # overall Runtime measuring
    dimensions = scfxn_ref15_space
    objective = design_with_config
    default_parameters = [val for k, val in ref15_weights]
    # "GP, GBRT, ET, RF"
    estimator = sys.argv[1]
    loss_value = sys.argv[2]

    # global save_and_exit

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
        acq_func_kwargs=acq_func_kwargs,
        # n_initial_points=300
    )

    # Warmstart Optimizer
    with open('wsres.pkl', 'rb') as h:
        wsres = pickle.load(h)
    
    y_0 = [x[loss_value] for x in wsres]
    x_0 = [x['weights'] for x in wsres]
    optimizer.tell(x_0, y_0)


    # MAIN optimiaztion Loop with custom result handling
    # using all available cores with Pool()
    optimizer_results = []

    result_buffer = {}
    jobs = []
    _DONE = False

    with get_context('spawn').Pool(processes=cpu_count(), maxtasksperchild=2) as tp:
        # maxtasksperchild for controlling memory usage
        BREAK_WAITING = False
        cached_config = None
        jobs_for_current_config = 0
        # For Checking outstanding jobs befor terminate.
        pending_jobs = cpu_count()

        def break_waiting():
            global BREAK_WAITING
            BREAK_WAITING = True

        def save_and_exit():
            global _DONE, result_buffer
            print("wait for all Processes to finish and close pool")
            print(len(active_children()))
            break_waiting()

            while not _DONE:
                for job in jobs:
                    if not job.ready():
                        time.sleep(60)
                        print(job._job, job._cache, job._event, job._pool)
                        job._event.set()  # tell process in same thread that is is done
                        # print(active_children())

                if all([job.ready() for job in jobs]):
                    print('ALL JOBS READY')
                    tp.close()
                    break

                if not active_children():
                    print('No More active_children')
                    _DONE = True
                    break
                time.sleep(5)

            print([job.ready() for job in jobs])
            print([(key, len(val)) for key, val in result_buffer.items()])
            print(len(result_buffer))
            print('SAVING')
            took = time.time() - start_time
            print("Took for ALL: {} to run"
                  .format(time.strftime("%H: %M: %S", time.gmtime(took))))
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
            _DONE = True

        my_configs = []
        def make_process():
            global calls
            global cached_config
            global jobs_for_current_config
            global pending_jobs

            # while not reached n_calls make either new config and first job for it
            # of instantiate job for existing config until @param runs_per_config.

            if cached_config and jobs_for_current_config < runs_per_config:
                config = cached_config
                jobs_for_current_config += 1
            else:  # make new config reset counter to 1
                config = optimizer.ask()
                optimizer.update_next()
                cached_config = config
                jobs_for_current_config = 1
            # instantiate job knowing which config it belongs to.
            my_configs.append(config)
            job = tp.apply_async(
                objective,
                (config,),
                callback=partial(_callback, config=config),
                error_callback=error_callback,
            )
            # Append map_result to jobs list
            jobs.append(job)
            calls += 1

        def _callback(map_res, config=None):
            """
            Whenever a process finishes it calls config_callback with its result
            and a new process can be run
            """
            global result_buffer, results, num_callbacks, calls, pending_jobs

            num_callbacks += 1
            print('CALLBACK: ', num_callbacks)
            # only frozenset and tuple are hashable
            # c_hash = hash(frozenset(config))
            # sorted(config)
            c_hash = str(sorted(config))

            # Check if key exists

            if c_hash in result_buffer.keys():

                if len(result_buffer[c_hash]) == (runs_per_config - 1):
                    # if this is the last run for config
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
                    # optimizer_results.append(opti_res)
                else:  # not last run yet, append to buffer for this config
                    print('APPEND map_res ')
                    result_buffer[c_hash].append(map_res)
            else: # otherwise make new buffer entry
                print('NEW BUFFER')
                result_buffer.update({c_hash: [map_res]})

		# while there are calls left on the budget
            if calls < n_calls:
                # Make New Process
                make_process()
            else:
                # wait for all processes to finish then exit Pool
                pending_jobs -= 1
                print('PENDING: ', pending_jobs)
                print('ACTIVE CHILDREN: ', len(active_children()))
                # print([job.ready() for job in jobs])

                # TODO: think HARD why there is one left job that is not ready ??
                # has not returned? counter error?
                if pending_jobs == 0:
                    save_and_exit()

        def error_callback(r):
            logger.log(level=logging.ERROR, msg=r)

        for _ in range(cpu_count()):
            # map runs_per_config processes with same config
            make_process()

        while not BREAK_WAITING:
            time.sleep(60)
            pass
