import logging

import pyrosetta as prs
import ray

from config import _init_method, _objective, dummy_objective


@ray.remote
class PRSActor(object):
    """docstring for PRSActor"""

    def __init__(self, arg):
        super(PRSActor, self).__init__()
        self.arg = arg
        prs.init()

    def __init__(self, initializer=None, initargs=None, idx=None):
        if initializer:
            initargs = initargs or ()
            _init_method(*initargs)
        self.idx = idx
        self.logger = logging.getLogger(__name__)

    def ping(self):
        # Used to wait for this actor to be initialized.
        # Also used to check if busy

        return self.idx

    def evaluate_config(self, config, run, pdb, target=None):
        self.logger.debug('Actor %d evaluates run %d', self.idx, run)

        return _objective(config, run, pdb, target=target)

    def run_batch(self, func, batch):
        results = []

        for args, kwargs in batch:
            args = args or ()
            kwargs = kwargs or {}

            if func is None:
                results.append(self.evaluate_config(*args, **kwargs))
            else:
                try:
                    results.append(func(*args, **kwargs))
                except Exception as e:
                    results.append(PoolTaskError(e))

        return results

    def shutdown(self):
        ray.actor.exit_actor()
