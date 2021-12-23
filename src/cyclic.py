
import argparse
import logging
import os
import sys
import time

import pandas as pd
import pyrosetta as prs
import ray

from bayesopt import BayesOpt
from distributor import Distributor
from manager import OptimizationManager

prs.init()
# make cluster recognize directory
sys.path.append(os.getcwd())


logger = logging.getLogger("APP")
logger.setLevel(logging.DEBUG)


if __name__ == "__main__":

    parser.add_argument(
        "-r_pw",
        "--redis_password",
        default=None,
        type=str,
        help="password to use for ray cluster redis authentication"
    )

    if args.redis_password:

        ray.init(address='auto', _redis_password=args.redis_password)
    else:
        ray.init(address='auto')
    print('''This cluster consists of
        {} nodes in total
        {} CPU resources in total
    '''.format(len(ray.nodes()), ray.cluster_resources()['CPU']))
    if args.cores == 0:  # no cores defined
        cores = ray.cluster_resources()['CPU']
    else:
        cores = args.cores
    logger.debug('Running on %d cores', cores)
    # # SIGNAL
    # signal = SignalActor.remote()
    pdb = args.pdb
    for i in range(cycles):
        # DISTRIBUTOR
        distributor = Distributor()
        # OPTIMIZER
        optimizer = BayesOpt()
        # MANAGER
        manager = OptimizationManager()
        manager.init(
            args.loss,
            pdb=pdb,
            distributor=distributor,
            optimizer=optimizer,
            estimator=args.estimator,
            identifier=args.id,
            test_run=args.test_run,
            n_cores=int(cores),
            evals=int(args.evals),
            rpc=int(args.runs_per_config),
            mtpc=int(args.max_tasks_per_child),
            cooldown=args.cooldown,
            out_dir=args.output_dir,
        )

        result = manager.run(report=True)
        winner_pose = prs.distributed.packed_pose.core.to_pose(
            result.groupby('run').mean().nsmallest(1, loss).pose)
        prs.dump_pdb(
            winner_pose, '../benchmark/allosteric/current_best'+str(i)+'_'+args.pdb+'.pdb')
        pdb = 'current_best'+str(i)+'_'+args.pdb+'.pdb'
