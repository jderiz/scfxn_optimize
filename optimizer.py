import logging
import os
import pickle
import random
import sys
import time

import numpy as np
import pandas as pd
from dask.distributed import as_completed
from dask_jobqueue import SLURMCluster
from skopt import Optimizer, Space, callbacks

import hyperparams
from relax import initialize, relax_with_config


class Optimizer:
    def __init__(
        self,
        loss,
        pdb=None,
        estimator="RF",  # "dummy" for random search
        identifier=None,
        test_run=False,
        n_calls=200,
        rpc=5,  # runs_per_config
        result_dir="results",
        warm_start=None,
        cores=None,
        mtpc=None,  # maxtasksperchild
        xi=0.01,  # starting value if cooldown
        kappa=1.69,  # starting value if cooldown
        cooldown=False,
        space_dimensions=None,
        save_pandas=True,
    ):
        # counter
        calls = 0
        # COOLDOWN LOOUP
        final_xi = 0.001
        final_kappa = 0.01
        _xi = xi
        _kappa = kappa
        xi_kappa_lookup = pd.DataFrame(None, index=range(n_calls))
        xi_kappa_lookup["iter"] = range(1, n_calls + 1)
        xi_kappa_lookup["geospace"] = np.geomspace(0.001, 1, num=n_calls)
        xi_kappa_lookup["xi"] = xi - xi_kappa_lookup.geospace * (xi - final_xi)
        xi_kappa_lookup["kappa"] = kappa - \
            xi_kappa_lookup.geospace * (kappa - final_kappa)
        
        # TEST CASE
        if test_run:
            print("DUMMY OBJECTIVE")
            objective = dummy_objective
            loss_value = 'ref15'
            init_method = None
            pandas = False
        else:
            objective = relax_with_config
            init_method = initialize


        # SEARCH SPACE

        if not space_dimensions:
            try:
                dimensions = Space.from_yaml("space.yml")
            except FileNotFoundError as e:
                print(
                    e,
                    "\n Could not read file space.yml, \n \
                        create it or supply another file path for searchspace creation",
                )
                pass
        else:
            with open(space_dimensions, "rb") as h:
                dimensions = pickle.load(h)
        # setup OPTIMIZER
        acq_func_kwargs = {"xi": xi, "kappa": kappa}
        optimizer = Optimizer(
            random_state=5,
            dimensions=dimensions,
            base_estimator=estimator,
            acq_func_kwargs=acq_func_kwargs,
            n_initial_points=rpc * 2,
        )

        # CLUSTER
        cluster = SLURMCluster()
        results = pd.DataFrame()

    def log_res(map_res: dict) -> None:
        # TODO: find type of future.result() in dask

        print(map_res)

        for res in map_res:
            res.update({"config": config})
            res.update({"c_hash": c_hash})
            res.update({"run": run})
        results.extend(map_res)
        optimizer.tell(config, (sum([x[self.loss]
                                     for x in map_res]) / len(map_res)))

    def run() -> None:
        # map initial runs
        starting_batches = self.cluster.map(
            partial(objective(pdb)), [config for _ in range(rpc)])
        seq = as_completed(starting_batches)

        for res in seq:
            self.log_res(res.result())

            if call <= n_calls:
                config = self.optimizer.ask()
                new_batch = self.cluster.submit(partial())
            seq.add(new_batch)
        pass

    def save_and_exit() -> bool:
        print(cluster)

        if self.pandas:
            df = pd.DataFrame(results)
            weights = df.config.apply(lambda x: pd.Series(x))
            weights.columns = ["fa_rep_" + str(num) for num in range(7)]
            results = pd.concat([df, weights], axis=1)
        with open(
            "results/{}_res_{}.pkl".format(identify, base_estimator),
            "wb",
        ) as file:
            pickle.dump(results, file)
        print("TERMINATING")

    def _cooldown() -> None:

        new_kappa = xi_kappa_lookup.kappa[self.calls]
        new_xi = xi_kappa_lookup.xi[self.calls]
        print("UPDATE XI KAPPA \n", new_xi, new_kappa)
        self.optimizer.acq_optimizer_kwargs.update(
            {"xi": new_xi, "kappa": new_kappa})
        self.optimizer.update_next()


# use for test purposes
def dummy_objective(pdb, config) -> dict:
    # print('TEST')
    # time.sleep(random.randint(5, 15))

    return {
        "bloss62": random.randint(1, 100),
        "ref15": random.randint(1, 50),
        "scfxn": random.randint(1, 46),
        "score": random.randint(1, 20),
    }
