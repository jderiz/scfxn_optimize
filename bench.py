import logging
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

logger = logging.getLogger("rosetta")

if __name__ == "__main__":
    """
    Main execution for optimization run.
    """
    # Setup result folder
    os.makedirs("results", exist_ok=True)

    # instantiate result array and specific number calls to objective per optimizer
    res = []
    n_calls = 50  # Objective Function evaluations
    start_time = time.time()  # overall Runtime measuring
    dimensions = scfxn_ref15_space
    # def test(x):
    #     time.sleep(5)
    #     print('\n \n \n \n \n')
    #     return x
    objective = design_with_config
    # objective = test
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

    # Run for n_calls asking #cores points each time

    run = 0

    # MAIN optimiaztion Loop with custom result handling
    # using all available cores with Pool()
    results_from_config = {}
    configs = {}
    semaphore_open = False
    calls = 0
    # results = []
    with get_context("spawn").Pool(processes=cpu_count()) as tp:
        print("initialized Worker Pool with max CPUS")

        def join_and_exit():
            print("wait for all Processes to finish and close pool")
            tp.join()
            tp.close()

        def make_config_batch():
            # while not reached n_calls

            if calls <= n_calls:
                config = optimizer.ask()
                # map 4 processes with same config calling callback once all are done
                r = tp.map_async(
                    objective,
                    [config for _ in range(runs_per_config)],
                    callback=partial(config_callback, config=config),
                    error_callback=error_callback,
                )
                r.wait()
                # results.append(r)
                calls += 1
            else:
                join_and_exit()

        def config_callback(res, config=None):
            """
            Whenever a process finishes it calls config_callback with its result
            and a new batch of processes can be run
            """
            print("map_async returning")
            _res = {
                "bloss62": sum([x["bloss62"] for x in res]) / len(res),
                "ref15": sum([x["ref15"] for x in res]) / len(res),
                "scfxn": sum([x["scfxn"] for x in res]) / len(res),
                "weights": configs[config],
            }
            print(_res)
            res.append(_res)
            optimizer.tell(config, _res["bloss62"])
            optimizer.update_next()  # check if necessary
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

        for _ in range(initial_batches):
            config = optimizer.ask()
            print(config)
            # results_from_config.update({hash(config): []})
            # map 4 processes with same config calling callback once all are done
            r = tp.map_async(
                objective,
                [config for _ in range(runs_per_config)],
                callback=partial(config_callback, config=config),
                error_callback=error_callback,
            )
            initial_map_results.append(r)

        for r in initial_map_results:
            r.wait()

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
