import argparse
from multiprocessing import Pool, cpu_count, get_context

import optimization
from relax import initialize

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
    parser.add_argument("-pdb",
                        type=str,
                        default=None,
                        help="ONLY DEV: specify which protein should be used")
    parser.add_argument(
        "-evals",
        type=int,
        default=200,
        help="number of points to evaluate")
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
        type=int,
        help="spcify number of cores to use if 0 then use all",
    )
    parser.add_argument(
        "-mtpc",
        "--max-tasks-per-child",
        default=3,
        type=int,
        help="limit the number of task a worker process can complete before respawning a new process ( useful for freeing RAM)",
    )
    parser.add_argument(
        "-range",
        default=0.25,
        type=float,
        help="how much the optimization can deviate from ref15 weights, defaults to 0.25")
    parser.add_argument(
        "-xi",
        default=0.01,
        type=float,
        help="optimizer argument to manage explore vs. exploit higher==> explore")
    parser.add_argument(
        "-kappa",
        default=1.69,
        type=float,
        help="optimier argument to manage explore vs. exploit higher==> explore")
    parser.add_argument(
        "-config",
        type=str,
        default=None,  # uses the protocolls default config
        help="If a config path to a pickled list, series or DataFrame that holds it is supplied this particular config is evaluated -evals times and the results are stored with all information.")
    parser.add_argument(
        "-id",
        default=None,
        type=str,
        help="Identity string for storing the results")
    parser.add_argument(
        "-no_struct",
        action="store_true",
        help="not saving the structures saves enormous amounts of space")
    parser.add_argument("-dict_out", action="store_true",
                        help="save result in dict format rather then pandas DataFrame")

    parser.add_argument("-cooldown", type=bool,
                        help='specifies if cooldown should be applied to the explore exploit params and if so linnear or logarithmic',
                        default=False, )
                        # choices=['lin', 'slow', 'fast'])
    args = parser.parse_args()
    print(args)

    if not args.cores:
        cores = cpu_count()
    else:
        cores = args.cores

    if args.config != None:
        #     # do design instead of optimization
        #     optimization.design(args.config, identify=args.id, evals=args.evals,
        #                         mtpc=args.max_tasks_per_child)
        # pass
        optimization.relax(
            identify=args.id, config_path=args.config, evals=args.evals, pdb=args.pdb)
    else:
        print('RUN Optimizer')
        optimization.init(
            args.loss,
            pdb=args.pdb,
            estimator=args.estimator,
            identifier=args.id,
            test_run=args.test_run,
            cores=int(cores),
            number_calls=int(args.evals),
            rpc=int(args.runs_per_config),
            mtpc=int(args.max_tasks_per_child),
            weight_range=args.range,
            xi=args.xi,
            kappa=args.kappa,
            cooldown=args.cooldown
        )
        optimization.start_optimization()
