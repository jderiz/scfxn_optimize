#!/usr/bin/env python3

import pickle
import pandas as pd
from multiprocessing import (Pipe, Pool, active_children, cpu_count,
                             current_process, get_context)

from joblib import Parallel, delayed
import pyrosetta as prs
from design import design_with_config

prs.init(
options="-ex1 -ex2", set_logging_handler=True, extra_options="-mute basic -mute core -mute protocols -linmem_ig 10 -archive_on_disk /tmp/rosetta"
)  # no output from the design process

if __name__ == "__main__":
    with open('all_sm.pkl', 'rb') as h:
        all_sm = pickle.load(h)

    config = all_sm.mean().values

    with get_context('spawn').Pool(processes=cpu_count(), maxtasksperchild=3) as tp:

        result_set = tp.map(design_with_config, [config for i in range(1000)] ) 
        with open('results/final_results.pkl') as h:
            res = pickle.load(h)
            res.append(pd.DataFrame(result_set))
            pickle.dump(res)
