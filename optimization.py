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
import numpy as np
import pandas as pd
from skopt import Optimizer, callbacks

import hyperparams
from design import design_with_config, initialize


def init(
    loss,
    estimator="RF",
    identifier=None,
    test_run=False,
    number_calls=200,
    rpc=5,  # runs_per_config
    result_dir="results",
    warm_start=None,
    cores=None,
    mtpc=None,  # maxtasksperchild
    weight_range=0.25,
    xi=0.01,  # starting value if cooldown
    kappa=1.69,  # starting value if cooldown
    cooldown=False,
    space_dimensions=None,
    save_pandas=True
):
    """
    Init Optimization 
    """

    # define global variables all functions can access
    global _DONE 
    global identify
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
    global runs_per_config 
    global loss_value 
    global base_estimator 
    global tp 
    global _cooldown 
    global pandas 
    global final_xi 
    global final_kappa 
    global _xi 
    global _kappa 
    global logger
    # Setup Logging
    logger = multiprocessing.log_to_stderr()
    logger.setLevel(logging.WARNING)
    # EXPLORE/EXPLOIT 
    final_xi = 0.001
    final_kappa = 0.01
    _xi = xi
    _kappa = kappa

    if not identifier:
        log_handler = logging.FileHandler('mp.log')
    else:
        log_handler = logging.FileHandler('mp_'+identifier+'.log')
    log_handler.setLevel(logging.DEBUG)
    logger.addHandler(log_handler)
    pandas = save_pandas

    identify = identifier
    # cached_config = None
    loss_value = loss
    # jobs_for_current_config = 0
    n_calls = number_calls
    runs_per_config = rpc
    calls = 0  # iteration counter
    num_callbacks = 0  # keep track of callback calls
    jobs = []
    _DONE = False
    result_buffer = {}
    results = []
    base_estimator = estimator
    _cooldown = cooldown

    # Setup result folder
    os.makedirs(result_dir, exist_ok=True)

    # Parallel handling

    if cores:
        _cores = cores
    else:
        _cores = cpu_count()
    # overall Runtime measuring
    start_time = time.time()

    # optimizer dimensions setup, either take default ref15 or pickled from user input

    if not space_dimensions:
        hyperparams.set_range(weight_range)
        dimensions = hyperparams.get_dimensions()
    else:
        with open(space_dimensions, 'rb') as h:
            dimensions = pickle.load(h)

    if test_run:
        print("DUMMY OBJECTIVE")
        objective = dummy_objective
        loss_value = 'ref15'
        init_method = None
        pandas = False
    else:
        objective = design_with_config
        init_method = initialize
    # Higher values for xi, kapp --> more exploration
    # DEFAULTS: xi:0.01, kappa:1.96
    # TODO: implement cooldown [start_values, end_values]
    acq_func_kwargs = {"xi": xi, "kappa": kappa}
    optimizer = Optimizer(
        random_state=5,
        dimensions=dimensions,
        base_estimator=base_estimator,
        acq_func_kwargs=acq_func_kwargs,
        # n_initial_points=300
    )
    # Warmstart Optimizer with previous results

    if warm_start:
        with open(warm_start, "rb") as h:
            wsres: pd.DataFrame = pickle.load(h)

        # we only want to tell the optimizer one y_0 per config
        c_group = wsres.goupby('c_hash')

        y_0 = [x[loss_value].mean() for x in c_group]
        x_0 = [x["weights"][0] for x in c_group]
        optimizer.tell(x_0, y_0)

    # always spawn:  https://pythonspeed.com/articles/python-multiprocessing/
    tp = get_context('spawn').Pool(
        cores, initializer=init_method, maxtasksperchild=mtpc)


def dummy_objective(config) -> dict:
    # print('TEST')
    # time.sleep(random.randint(5, 15))

    return {
        "bloss62": random.randint(1, 100),
        "ref15": random.randint(1, 50),
        "scfxn": random.randint(1, 46),
    }


def save_and_exit() -> None:
    """All Callbacks have returned, tell all remaining jobs they are done, terminate pool and save results"""
    print(len(active_children()))
    global jobs
    global _DONE
    global base_estimator
    global tp
    global identify

    # tell pool its about to terminate
    tp.close()
    print(tp._state)
    print([job.ready() for job in jobs])

    print("SAVING")
    took = time.time() - start_time
    print(
        "Took for ALL: {} to run".format(
            time.strftime("%H: %M: %S", time.gmtime(took))
        )
    )
    # save custom results, as pandas or dict

    if pandas:
        df = pd.DataFrame(results)
        weights = df.config.apply(lambda x: pd.Series(x))
        weights.columns = [name for name, _ in hyperparams.ref15_weights]
        results = pd.concat([df, weights], axis=1)
    with open(
        "results/{}_{}_res_{}.pkl".format(
            time.strftime("%H_%M"), base_estimator, identify
        ),
        "wb",
    ) as file:
        pickle.dump(results, file)
    print("TERMINATING")
    tp.terminate()


def _callback(map_res, config=None, run=None) -> None:
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
    global _DONE
    num_callbacks += 1
    # print(map_res)
    # add weights to each entry as well as config hash for later analysis
    print("CALLBACK: ", num_callbacks)
    c_hash = str(sorted(config))

    for res in map_res:
        res.update({'config': config})
        res.update({'c_hash': c_hash})
        res.update({'run': run})
    # print(map_res)
    results.extend(map_res)
    optimizer.tell(
        config, (sum([x[loss_value] for x in map_res]) / len(map_res))
    )

    # make n_calls batches

    if calls < n_calls:
        print(calls, n_calls)
        # Make New Process batch
        make_batch()
    else:
        if len(results) == n_calls*runs_per_config:
            print('got ', len(results), ' results from ', len(jobs), ' jobs from which ', len(
                [job for job in jobs if job.ready()]), 'report ready')
            _DONE = True  # release waiting Loop
            save_and_exit()


def error_callback(r):
    global logger
    logger.log(level=logging.ERROR, msg=r)


def make_batch(config=None) -> None:
    # global cached_config
    # global jobs_for_current_config
    global runs_per_config
    global calls
    global tp
    # global n_calls

    calls += 1
    print('Make Batch: ', calls)
    
    # scale xi & kapp(explore/exploite) according to chosen curve
    if _cooldown:
        # either lin, log(fast), neg exp(slow)
        # TODO: contract/expand cooldown
        # TODO: precompute and store in optimizer singleton, make dependent on start /end values
        global final_kappa, final_xi, n_calls
        if _cooldown == 'lin':
            new_xi = _xi - (((_xi-final_xi)/n_calls)*calls)
            new_kappa = _kappa - (((_kappa-final_kappa)/n_calls)*calls)
        elif _cooldown == 'fast':
            new_xi = np.logspace(0.1, -3, num=n_calls, endpoint=True)[calls]
            new_kappa = np.logspace(1, -1, num=n_calls, endpoint=True)[calls]
        elif _cooldown == 'slow':
            # TODO: make proper
            df = pd.DataFrame(None, index=range(n_calls))
            df[1] = range(1, n_calls+1)
            df['xi'] = df[1].apply(lambda x: _xi - (x/n_calls*(x*-(final_xi-_xi)))/n_calls )
            df['kappa'] = df[1].apply(lambda x: _kappa - (x/n_calls*(x*-(final_kappa-_kappa)))/n_calls )
            new_xi = df.xi[calls]
            new_kappa = df.kappa[calls]


        print('UPDATE XI KAPPA \n', new_xi, new_kappa)
        optimizer.acq_optimizer_kwargs.update(
            {'xi': new_xi, 'kappa': new_kappa})
        optimizer.update_next()

    if not config:
        config = optimizer.ask()
    job = tp.map_async(
        objective,
        [config for _ in range(runs_per_config)],
        callback=partial(_callback, config=config, run=calls),
        error_callback=error_callback,
    )
    # Append map_result to jobs list
    jobs.append(job)
    # increase calls counter


def design(config_path=None, identify=None, evals=1000, mtpc=3, cores=cpu_count()):
    """
        Do an actual design run with a single config.
        The config either needs to be a list with weight values in 
        "correct" order. Or a pd.Series or DataFrame object with corresponding 
        column names.
    """
    print(config_path)

    if config_path == 'ref15':
        # use ref15 as default weights
        config = 'ref15'  # [weight for _, weight in hyperparams.ref15_weights]
    else:
        with open(config_path, 'rb') as h:
            config = pickle.load(h)

            if type(config) == list:
                pass
            elif (type(config) == pd.DataFrame) and (len(config) == 1):
                c = []

                for w in [name for name, _ in hyperparams.ref15_weights]:
                    c.append(config.iloc[0][w])
                config = c
                print(config)
            elif type(config) == pd.Series:
                c = []

                for w in [name for name, _ in hyperparams.ref15_weights]:
                    c.append(config[w])
                    # TODO: handle series access by label.
                config = c
            else:
                exit(TypeError(
                    'the config must be either in list or DataFrame (1, len(weights)) format'))

    with get_context('spawn').Pool(processes=cores, initializer=initialize, maxtasksperchild=mtpc) as tp:
        result_set = tp.map(design_with_config, [config for i in range(evals)])
        # TODO: refactor into helper function
        res = pd.DataFrame(result_set)
        with open('results/design_{}_{}.pkl'.format(evals, identify), 'wb') as h:
            pickle.dump(res, h)


def start_optimization():
    """
    start the optimization process, and wait for it to finish
    """
    global _DONE
    global calls
    global runs_per_config
    _DONE = False
    ref15_config = [val for k, val in hyperparams.get_ref15()]
    # initial runs, starting with ref15

    make_batch(ref15_config)

    for _ in range(int(_cores/runs_per_config) - 1):
        make_batch()

    while not _DONE:
        # print('Wait...')
        time.sleep(5)
        pass
