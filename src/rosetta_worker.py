import logging

import pyrosetta as prs
import ray
from ray.util import inspect_serializability

# import objective as it cannot be pickled
from config import _init_method, _objective


@ray.remote
class PRSActor(object):
    """docstring for PRSActor"""

    def __init__(self, arg):
        super(PRSActor, self).__init__()
        self.arg = arg
        prs.init()

    def __init__(self, initargs=None, idx=None):
        """ Initialize the actor """

        if _init_method:  # if initis set, call it with initargs
            initargs = initargs or ()
            _init_method(*initargs)
        self.idx = idx
        self.logger = logging.getLogger(__name__)
        # self.logger.debug('CHECK SERIALIZABLE OBJECTIVES FUNCTIONS')
        # inspect_serializability(_init_method, 'Initializer Method')
        # inspect_serializability(_objective, 'Objective Method')

    def ping(self):
        """
        Ping an actor instance.
        Used to check if the actor is still alive.
        Also used to check if actor was initialized and is ready to receive work
        """

        return self.idx

    def evaluate_config(self, config, fargs, run) -> dict:
        """
        evaluates configuration and returns a list with the score and the
        run this evaluation belongs to.
        """
        self.logger.debug('Actor %d evaluates run %d for %s',
                          self.idx, run, fargs)
        res: dict = _objective(config, *fargs)
        res.update({'run': run})

        return res

    def shutdown(self):
        """
        shutdown the actor
        """
        ray.actor.exit_actor()
