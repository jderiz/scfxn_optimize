import logging
import multiprocessing
import os
import pickle
import random
import sys
import time
from functools import partial
from multiprocessing import (Pipe, Pool, active_children, cpu_count,
                             current_process, get_context)

from numpy import repeat

from design import design_with_config, initialize
from hyperparams import ref15_weights, scfxn_ref15_space
from skopt import (Optimizer, callbacks, forest_minimize, gbrt_minimize,
                   gp_minimize)

# Setup Logging
logger = multiprocessing.log_to_stderr()
logger.setLevel(logging.WARNING)


def init(
    loss,
    base_estimator="RF",
    test_run=False,
    n_calls=200,
    runs_per_config=5,
    result_dir="results",
    warm_start=None,
    cores=None,
):
    """
    Init Optimization 
    """
    # define global variables all functions can access
    global _DONE
    global _cores
    global result_buffer
    global results
    global n_calls
    global calls
    global num_callbacks
    global jobs
    global start_time
    global objective
    global optimizer
    global cached_config
    global runs_per_config
    global loss_value
    global jobs_for_current_config
    global base_estimator

    cached_config = None
    loss_value = loss
    jobs_for_current_config = 0
    n_calls = n_calls
    runs_per_config = runs_per_config
    calls = 0  # iteration counter
    num_callbacks = 0  # keep track of callback calls
    jobs = []
    _DONE = False
    result_buffer = {}
    results = []
    base_estimator = base_estimator

    # Setup result folder
    os.makedirs(result_dir, exist_ok=True)

    # Parallel handling

    if cores:
        _cores = cores
    else:
        _cores = cpu_count()
    # overall Runtime measuring
    start_time = time.time()

    # Optimizer setup
    dimensions = scfxn_ref15_space

    if test_run:
        print("DUMMY OBJECTIVE")
        objective = dummy_objective
    else:
        objective = design_with_config
    # more exploit then explore: DEFAULTS: xi:0.01, kappa:1.96
    # TODO: implement cooldown [start_values, end_values]
    xi = 0.001
    kappa = 0.1
    acq_func_kwargs = {"xi": xi, "kappa": kappa}
    optimizer = Optimizer(
        dimensions=dimensions,
        base_estimator=base_estimator,
        acq_func_kwargs=acq_func_kwargs,
        # n_initial_points=300
    )
    # Warmstart Optimizer

    if warm_start:
        with open(warm_start + ".pkl", "rb") as h:
            wsres = pickle.load(h)

        y_0 = [x[loss_value] for x in wsres]
        x_0 = [x["weights"] for x in wsres]
        optimizer.tell(x_0, y_0)


def dummy_objective(config):
    # print('TEST')
    time.sleep(1)

    return {
        "bloss62": random.randint(1, 100),
        "ref15": random.randint(1, 50),
        "scfxn": random.randint(1, 46),
    }


def save_and_exit(tp):
    print("wait for all Processes to finish and close pool")
    print(len(active_children()))
    global jobs
    global _DONE
    global base_estimator

    while not _DONE:
        for job in jobs:
            if not job.ready():
                time.sleep(60)
                print(job._job, job._cache, job._event, job._pool)
                job._event.set()  # tell process in same thread that is is done
                # print(active_children())

        if all([job.ready() for job in jobs]):
            print("ALL JOBS READY")
            tp.close()

            break

        if not active_children():
            print("No More active_children")
            _DONE = True

            break
        time.sleep(5)

    print([job.ready() for job in jobs])
    print([(key, len(val)) for key, val in result_buffer.items()])
    print(len(result_buffer))
    print("SAVING")
    took = time.time() - start_time
    print(
        "Took for ALL: {} to run".format(
            time.strftime("%D_%H: %M: %S", time.gmtime(took))
        )
    )
    # save custom results
    with open(
        "results/{}_{}_res_{}_{}.pkl".format(
            time.strftime("%D_%H_%M"), base_estimator, str(n_calls), loss_value,
        ),
        "wb",
    ) as file:
        pickle.dump(results, file)
    print("TERMINATING")
    tp.terminate()
    _DONE = True


def make_process(tp, config=None):
    # TODO: lock

    # while not reached n_calls make either new config and first job for it
    # of instantiate job for existing config until @param runs_per_config.
    global cached_config
    global jobs_for_current_config
    global runs_per_config
    global calls

    if config:
        config = config
        cached_config = config
        jobs_for_current_config += 1
    elif cached_config and jobs_for_current_config < runs_per_config:
        config = cached_config
        jobs_for_current_config += 1
    else:  # make new config reset counter to 1
        config = optimizer.ask()
        optimizer.update_next()
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

    # TODO: release lock


def _callback(map_res, config=None):
    """
    Whenever a process finishes it calls config_callback with its result
    and a new process can be run. 
    """
    # TODO: lock
    global num_callbacks
    global loss_value
    global result_buffer
    global results
    

    num_callbacks += 1
    print("CALLBACK: ", num_callbacks)
    # only frozenset and tuple are hashable
    # c_hash = hash(frozenset(config))
    # sorted(config)
    c_hash = str(sorted(config))

    # Check if key exists

    if c_hash in result_buffer.keys():

        if len(result_buffer[c_hash]) == (runs_per_config - 1):
            # if this is the last run for config
            print("SAVE RESULT")
            # if so, compute mean over results for config and tell optimizer
            result_buffer[c_hash].append(map_res)
            # print(result_buffer[c_hash])
            res = result_buffer[c_hash]
            _res = {
                "sequence": map_res["sequence"],
                "bloss62": sum([x["bloss62"] for x in res]) / len(res),
                "ref15": sum([x["ref15"] for x in res]) / len(res),
                "scfxn": sum([x["scfxn"] for x in res]) / len(res),
                "weights": config,
            }
            results.append(_res)
            optimizer.tell(config, _res[loss_value])
            # del result_buffer[c_hash]  # dont need that buffer
            # optimizer_results.append(opti_res)
        else:  # not last run yet, append to buffer for this config
            print("APPEND map_res ")
            result_buffer[c_hash].append(map_res)
    else:  # otherwise make new buffer entry
        print("NEW BUFFER")
        result_buffer.update({c_hash: [map_res]})

    # while there are calls left on the budget

    # TODO: release lock

    if calls < n_calls:
        # Make New Process
        make_process()
    else:
        # wait for all processors to finish then exit Pool

        while active_children():
            print(len(jobs))
            time.sleep(60)

        # pending_jobs -= 1
        # print("PENDING: ", pending_jobs)
        # print("ACTIVE CHILDREN: ", len(active_children()))
        # # print([job.ready() for job in jobs])

        # # TODO: think HARD why there is one left job that is not ready ??
        # # has not returned? counter error?

        # if pending_jobs == 0:
        #     save_and_exit()


def error_callback(r):
    logger.log(level=logging.ERROR, msg=r)


def start():
    """
    start the optimization process, and wait for it to finish
    """
    _DONE = False
    ref15_config = [val for k, val in ref15_weights]
    # initial runs

    make_process(ref15_config)

    for _ in range(_cores):
        make_process()

    while not _DONE:
        time.sleep(60)
        pass
