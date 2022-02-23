import argparse
import logging
import os
import sys
import time

import pyrosetta as prs
import ray

from bayesopt import BayesOpt
from config import result_path
from distributor import Distributor
from manager import OptimizationManager

# make cluster recognize directory
# sys.path.append(os.getcwd())

logger = logging.getLogger("APP")
logger.setLevel(logging.DEBUG)


if __name__ == "__main__":
    """

    Main execution script for Optimization run

    """

    ##############################
    # ARG PARSER
    ###############################
    description = "Run optimization of Energy Function weights for some objective\n\
                    function in parallel."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "-r_pw",
        "--redis_password",
        default=None,
        type=str,
        help="password to use for ray cluster redis authentication"
    )
    parser.add_argument(
        "-e",
        "--estimator",
        default="RF",
        type=str,
        help="the base_estimator to be used by the Optimizer in skopt (RF, ET, GBRT, GP) \
        or dummy for randomsearch",
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
                        help="specify which protein should be used")
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
        type=str,
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
        help="If a config path to a pickled list, series or DataFrame that holds it is supplied this particular config is evaluated -evals times and the results are stored with all information.")
    parser.add_argument(
        "-id",
        default=None,
        type=str,
        help="Identity string for storing the results")
    parser.add_argument(
        "-target",
        default=None,
        type=str,
        help='Pdb file to compare against'
    )
    parser.add_argument(
        "-cycles",
        default=1,
        type=int,
        help="how many optimization cycles, each new cycle uses best observed so far as starting pose"
    )
    parser.add_argument(
        "-no_struct",
        action="store_true",
        help="not saving the structures saves enormous amounts of space")
    parser.add_argument("-dict_out", action="store_true",
                        help="save result in dict format rather then pandas DataFrame")

    parser.add_argument(
        "-cooldown",
        action='store_true',
        help='specifies if cooldown should be applied to the explore exploit params and if so linnear or logarithmic',
        default=False
    )
    parser.add_argument(
        "-crb",
        "--complete-run-batch",
        default=False,
        action="store_true")
    args = parser.parse_args()
    print(args)

    ########################
    # SETUP RAY CLUSTER
    ########################

    if args.redis_password:  # On HPC we need redis_password for inter-node comm

        ray.init(address='auto', _redis_password=args.redis_password)
    else:
        ray.init(address='auto')
    print('''This cluster consists of
        {} nodes in total
        {} CPU resources in total
    '''.format(len(ray.nodes()), ray.cluster_resources()['CPU']))
    if args.cores == 0:  # no cores defined,
        cores = ray.cluster_resources()['CPU']
    else:
        cores = args.cores

    logger.debug('Running on %d cores', cores)

    pdb = args.pdb
    target = args.target
    prot_name = pdb.split("_")[-1]  # last element is pdb name with ending
    #########################
    # SETUP Modules
    #########################
    # DISTRIBUTOR
    distributor = Distributor()
    # OPTIMIZER
    optimizer = BayesOpt()
    # MANAGER
    manager = OptimizationManager()
    manager.init(
        args.loss,
        pdb=pdb,
        target=target,
        distributor=distributor,
        optimizer=optimizer,
        estimator=args.estimator,
        identifier=args.id,
        test_run=args.test_run,
        n_cores=int(cores),
        evals=int(args.evals),
        rpc=int(args.runs_per_config),
        cooldown=args.cooldown,
        out_dir=args.output_dir,
    )

    ###################################
    # Optimization Cycle(s)
    ###################################
    for cycle in range(args.cycles):
        if args.config != None:  # when a specific config supplied dont optimize
            manager.no_optimize(
                identify=args.id,
                config_path=args.config,
                evals=args.evals,
                pdb=args.pdb)
            break
        elif args.cycles == 1:  # when only once cycle
            manager.run(complete_run_batch=args.complete_run_batch)
            break  # make sure we leave the loop
        else:
            manager.set_pdb(pdb)
            manager.set_cycle(cycle=cycle)
            # reinitialize optimizer, as we need a new one.
            manager.init_optimizer()
            result = manager.run(
                report=True, complete_run_batch=args.complete_run_batch)

            # choose absolute best, group mean is important for
            # optimizer to have some certainty of config quality
            # and not just a lucky one
            winner_pose = prs.distributed.packed_pose.core.to_pose(
                result.nsmallest(1, args.loss).pose.iloc[0])
            # write current_best to disk
            prs.dump_pdb(
                winner_pose, result_path+'current_best_cycle_'+str(i)+'_'+prot_name)
            # make pdb_path string for next iteration point to current best
            pdb = result_path+'current_best_cycle_'+str(i)+'_'+prot_name
            logger.info('FINISHED RUN %d saving result  at %s', i, pdb)

    #############################
    # SAVE AND EXIT
    #############################

    # We are done so save the result and terminate the cluster/Pool
    manager.save()
    distributor.terminate_pool()
    logger.debug('FINISHED OPTIMIZATION')
