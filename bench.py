import os
import pickle
import sys
import time

from joblib import Parallel, delayed
from skopt import (Optimizer, callbacks, forest_minimize, gbrt_minimize,
                   gp_minimize)

from design import design_with_config
from hyperparams import ref15_weights, scfxn_ref15_space

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
    objective = design_with_config
    default_parameters = [val for k, val in ref15_weights]
    # setup callbacks for logging
    timer_callback = callbacks.TimerCallback()
    forest_check = callbacks.CheckpointSaver(".forest_checkpoints.gz")
    gbrt_check = callbacks.CheckpointSaver(".gbrt_checkpoints.gz")
    gp_check = callbacks.CheckpointSaver(".gp_checkpoints.gz")
    # "GP, GBRT, ET, RF"
    estimator = sys.argv[1]
    cores = 64
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

    for i in range(n_calls):
        x = optimizer.ask(n_points=cores)
        y = Parallel(n_jobs=cores)(delayed(objective)(v) for v in x)
        optimizer.tell(x, y)
    res = optimizer.get_result()
    took = time.time() - start_time
    print("Took: {} to run".format(time.strftime("%H: %M: %S", time.gmtime(took))))
    with open("results/{}_res_{}.pkl".format(estimator, str(n_calls)), "wb") as file:
        pickle.dump(res, file)
