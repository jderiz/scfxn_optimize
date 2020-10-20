import pickle
import sys
import time
import os
from skopt import forest_minimize, gp_minimize, gbrt_minimize, callbacks

from skopt.plots import (plot_convergence, plot_evaluations, plot_objective)

from design import design_with_config
from hyperparams import scfxn_ref15_space, ref15_weights

if __name__ == "__main__":
    """
    Main execution for optimization run.
    """
    # Setup result folder
    os.makedirs('results', exist_ok=True)

    # instantiate result array and specific number calls to objective per optimizer
    res = []
    n_calls = 100  # Objective Function evaluations
    start_time = time.time()  # overall Runtime measuring

    dimensions = scfxn_ref15_space

    objective = design_with_config
    default_parameters = [val for k, val in ref15_weights]
    # setup callbacks for logging
    timer_callback = callbacks.TimerCallback()
    forest_check = callbacks.CheckpointSaver('.forest_checkpoints.gz')
    gbrt_check = callbacks.CheckpointSaver('.gbrt_checkpoints.gz')
    gp_check = callbacks.CheckpointSaver('.gp_checkpoints.gz')
    optimizer = sys.argv[1]
    xi = 0.001
    kappa = 0.1
    print(
        '_________start optimize________' +
        '_____________{}________________'.format(optimizer)
    )
    if optimizer in ['forest', 'all']:
        print('FOREST')
        forest_result = forest_minimize(
            func=objective,
            dimensions=dimensions,
            acq_func='PI',  # expected Improvement (evaluate whats best to use)
            n_calls=n_calls,
            x0=default_parameters,
            n_jobs=-1,
            xi=xi,
            kappa=kappa,
            callback=[timer_callback, forest_check]
        )  # All availiable Cores if aquisition =lbfgs
        with open("results/res_{}_{}_similar.pkl".format('forest', str(n_calls)), "wb") as file:
            pickle.dump(forest_result, file)
    if optimizer in ['gbrt', 'all']:
        print('GradientBoostedTrees')
        gbrt_res = gbrt_minimize(
            func=objective,
            dimensions=dimensions,
            acq_func='PI',
            n_calls=n_calls,
            x0=default_parameters,
            n_jobs=-1,
            xi=xi,
            kappa=kappa,
            callback=[timer_callback, gbrt_check]
        )
        with open("results/res_{}_{}_similar.pkl".format('gbrt', str(n_calls)), "wb") as file:
            pickle.dump(gbrt_res, file)
    if optimizer in ['gp', 'all']:
        print('Gaussian Processes')
        gp_res = gp_minimize(
            func=objective,
            dimensions=dimensions,
            acq_func="PI",
            n_calls=n_calls,
            n_jobs=-1,
            x0=default_parameters,
            xi=xi,
            kappa=kappa,
            callback=[timer_callback, gp_check]
        )
        with open("results/res_{}_{}_similar.pkl".format('gp', str(n_calls)), "wb") as file:
            pickle.dump(gp_res, file)

    # plot_convergence(*res)
    # plt.show()
    # plt.savefig('results/convergence.png', dpi=300, bbox_inches="tight")

    took = time.time() - start_time
    print("Took: {} to run".format(time.strftime("%H: %M: %S", time.gmtime(took))))
