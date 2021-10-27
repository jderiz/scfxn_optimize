import pyrosetta as prs
import ray
from config import _objective, _init_method

@ray.remote
class PRSActor(object):
    """docstring for PRSActor"""
    def __init__(self, arg):
        super(PRSActor, self).__init__()
        self.arg = arg
        self.ref_scfxn = get_fa_scorefxn()
        initialize()

    def evaluate_config(self, config, run, pdb):
        objective(config, run, pdb)

        
    def __init__(self, initializer=None, initargs=None):
        if initializer:
            initargs = initargs or ()
            initializer(*initargs)

    def ping(self):
        # Used to wait for this actor to be initialized.
        pass

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

