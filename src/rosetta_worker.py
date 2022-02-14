import logging

import pyrosetta as prs
import ray

# import objective as it cannot be pickled
from config import _init_method, _objective


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

    def shutdown(self):
        ray.actor.exit_actor()
