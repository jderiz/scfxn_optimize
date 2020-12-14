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
    estimator="RF",
    test_run=False,
    number_calls=200,
    rpc=5,
    result_dir="results",
    warm_start=None,
    cores=None,
    mtpc=None,
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
    global tp

    cached_config = None
    loss_value = loss
    jobs_for_current_config = 0
    n_calls = number_calls
    runs_per_config = rpc
    calls = 0  # iteration counter
    num_callbacks = 0  # keep track of callback calls
    jobs = []
    _DONE = False
    result_buffer = {}
    results = []
    base_estimator = estimator

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
        init_method = None
        n_calls = 200
    else:
        objective = design_with_config
        init_method = initialize
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

        # TODO: group by c_hash and then mean() over loss for each group
        y_0 = [x[loss_value] for x in wsres]
        x_0 = [x["weights"] for x in wsres]
        optimizer.tell(x_0, y_0)

    tp = Pool(cores, initializer=init_method, maxtasksperchild=mtpc)


def dummy_objective(config) -> dict:
    # print('TEST')
    time.sleep(random.randint(1, 5))

    return {
        "bloss62": random.randint(1, 100),
        "ref15": random.randint(1, 50),
        "scfxn": random.randint(1, 46),
    }


def save_and_exit() -> None:
    print("wait for all Processes to finish and close pool")
    print(len(active_children()))
    global jobs
    global _DONE
    global base_estimator
    global tp

    while not _DONE:
        for job in jobs:
            if not job.ready():
                time.sleep(5)
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
    print("SAVING")
    took = time.time() - start_time
    print(
        "Took for ALL: {} to run".format(
            time.strftime("%H: %M: %S", time.gmtime(took))
        )
    )
    # save custom results
    with open(
        "results/{}_{}_res_{}_{}.pkl".format(
            time.strftime("%H_%M"), base_estimator, str(n_calls), loss_value,
        ),
        "wb",
    ) as file:
        pickle.dump(results, file)
    print("TERMINATING")
    tp.terminate()
    _DONE = True


def _callback(map_res, config=None) -> None:
    """
    Whenever a process finishes it calls _callback with its result
    and a new process can be run. 
    """
    global num_callbacks
    global loss_value
    global results
    global tp
    global n_calls
    global calls
    num_callbacks += 1
    print("CALLBACK: ", num_callbacks)
    # print(map_res)
    # add weights to each entry
    c_hash = str(sorted(config))
    for res in map_res:
        res.update({'weights': config})
        res.update({'c_hash': c_hash})
    # print(map_res)
    results.extend(map_res)
    optimizer.tell(
        config, (sum([x[loss_value] for x in map_res]) / len(map_res))
    )

    if calls <= n_calls:
        print(calls, n_calls)
        # Make New Process batch
        make_batch()
    else:
        if num_callbacks == n_calls:
            save_and_exit()


def error_callback(r):
    logger.log(level=logging.ERROR, msg=r)


def make_batch(config=None) -> None:
    print('Make Batch')

    global cached_config
    global jobs_for_current_config
    global runs_per_config
    global calls
    global tp
    global n_calls

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
    # instantiate jobs knowing which config they belongs to.
    job = tp.map_async(
        objective,
        [config for _ in range(runs_per_config)],
        callback=partial(_callback, config=config),
        error_callback=error_callback,
    )
    # Append map_result to jobs list
    jobs.append(job)
    calls += 1


def start():
    """
    start the optimization process, and wait for it to finish
    """
    global _DONE
    global calls
    global runs_per_config
    _DONE = False
    ref15_config = [val for k, val in ref15_weights]
    # initial runs

    make_batch(ref15_config)
    calls += 1

    for _ in range(int(_cores/runs_per_config) - 1):
        make_batch()

    while not _DONE:
        print('Wait...')
        time.sleep(25)
        pass