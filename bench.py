import pickle
import time
import os
import matplotlib.pyplot as plt

from skopt import forest_minimize, gp_minimize, gbrt_minimize, callbacks
from skopt.plots import (  # only needed for bench black-box functions
    plot_convergence, plot_evaluations, plot_objective)

from design import Designer, design_with_config
from hyperparams import scfxn_ref15_space

if __name__ == "__main__":
    """
    Main execution for optimization run.
    """
    # Setup result folder
    os.makedirs('results', exist_ok=True)

    # instantiate result array and specific number calls to objective per optimizer
    res = []
    n_calls = 20  # Objective Function evaluations
    start_time = time.time()  # overall Runtime measuring

    dimensions = scfxn_ref15_space

    designer = Designer(filename='benchmark/1K9P_A_relax_0001.pdb')
    objective = design_with_config

    # setup callbacks for logging
    timer_callback = callbacks.TimerCallback()
    forest_check = callbacks.CheckpointSaver('.forest_checkpoints.gz')
    gbrt_check = callbacks.CheckpointSaver('.gbrt_checkpoints.gz')
    gp_check = callbacks.CheckpointSaver('.gp_checkpoints.gz')

    print(
        '_________start optimize________'
    )
    print('FOREST')
    forest_result = forest_minimize(
        func=objective,
        dimensions=dimensions,
        acq_func='EI',  # expected Improvement (evaluate whats best to use)
        n_calls=n_calls,
        # x0=default_parameters,
        n_jobs=-1,
        callback=[timer_callback, forest_check]
    )  # All availiable Cores if aquisition =lbfgs
    res.append(("forest", forest_result))

    print('GradientBoostedTrees')
    gbrt_res = gbrt_minimize(
        func=objective,
        dimensions=dimensions,
        acq_func='EI',
        n_calls=n_calls,
        n_jobs=-1,
        callback=[timer_callback, gbrt_check]
    )

    res.append(('gbrt', gbrt_res))
    print('Gaussian Processes')
    gp_res = gp_minimize(
        func=objective,
        dimensions=dimensions,
        acq_func="EI",
        n_calls=n_calls,
        n_jobs=-1,
        callback=[timer_callback, gp_check]
    )
    res.append(('gp', gp_res))

    with open("results/res_{}.pkl".format(str(n_calls)), "wb") as file:
        pickle.dump(res, file)
    plot_convergence(*res)
    plt.show()
    plt.savefig('results/convergence.png', dpi=300, bbox_inches="tight")

    took = time.time() - start_time
    print("Took: {} to run".format(time.strftime("%H: %M: %S", time.gmtime(took))))
