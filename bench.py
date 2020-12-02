import argparse
from multiprocessing import Pool, cpu_count

from optimization import Optimization

if __name__ == "__main__":
    """

    Main execution script for Optimization run

    """
    description = "Run optimization of Energy Function weights for some objective "
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-e",
        "--estimator",
        default="RF",
        help="the base_estimator to be used by the Optimizer",
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
        help="specifies a pickled object that holds information of previous runs and can be used to warm start the optimization by telling the optimiter about those points",
    )
    parser.add_argument(
        "-c",
        "--cores",
        default=0,
        help="spcify number of cores to use if 0 then use all",
    )

    args = parser.parse_args()
    print(args)

    if not args.cores:
        cores = cpu_count()
    else:
        cores = args.cores

    with Pool(processes=cores, maxtasksperchild=5) as tp:  # get_context('spawn')
        optimization = Optimization(
            tp,
            args.loss,
            base_estimator=args.estimator,
            test_run=args.test_run,
        )
        optimization.start()
