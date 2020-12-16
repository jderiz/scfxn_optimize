import argparse
from multiprocessing import Pool, cpu_count, get_context

import optimization
from design import initialize

if __name__ == "__main__":
    """

    Main execution script for Optimization run

    """
    description = "Run optimization of Energy Function weights for some objective\n\
                    function in parallel."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-e",
        "--estimator",
        default="RF",
        help="the base_estimator to be used by the Optimizer in skopt (RF, ET, GBRT, GP)",
    )
    parser.add_argument(
        "-l",
        "--loss",
        help="the dictionary entry resturned by the objective to use for optimization",
    )
    parser.add_argument(
        "-t",
        "--test-run",
        action="store_true",
        help="can be used as switch to evaluate dummy objctive",
    )
    parser.add_argument("-evals", default=200, help="number of points to evaluate")
    parser.add_argument(
        "-rpc",
        "--runs-per-config",
        default=5,
        help="how many evaluations per configuration",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="results",
        help="directory where the results get saved",
    )
    parser.add_argument(
        "-w",
        "--warm-start-file",
        help="specifies a pickled object that holds information of previous runs and can be used to warm start the optimization by telling the optimizer about those points",
    )
    parser.add_argument(
        "-c",
        "--cores",
        default=0,
        help="spcify number of cores to use if 0 then use all",
    )
    parser.add_argument(
        "-mtpc",
        "--max-tasks-per-child",
        default=3,
        help="limit the number of task a worker process can complete before respawning a new process ( useful for freeing RAM)",
    )

    args = parser.parse_args()
    print(args)

    if not args.cores:
        cores = cpu_count()
    else:
        cores = args.cores

    optimization.init(
        args.loss,
        estimator=args.estimator,
        test_run=args.test_run,
        cores=int(cores),
        number_calls=int(args.evals),
        rpc=int(args.runs_per_config),
        mtpc=int(args.max_tasks_per_child)
    )
    optimization.start()
