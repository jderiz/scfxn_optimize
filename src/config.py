"""
CONFIGURATION
Import your objective function. If you have an initializer function to run on 
worker instantiation also import it here and assign it to the corresponding
variables.
Set custom for xi and kappa if you whish else leave them at the default.
Specify a path where your results should get saved to.

"""
import logging

from allosteric import initialize, relax_with_config

# LOG LEVEL
level = logging.DEBUG
handler = logging.StreamHandler()
logging.getLogger().addHandler(handler)
logging.getLogger().setLevel(level)
# OBJECTIVE FUNCTION TO EVALUATE
_objective = relax_with_config
# SETUP THAT IS BEING DONE ONLY ONCE ON EACH CPU e.g initializing Rosetta, loading large files
_init_method = initialize
# SPACE DIMENSIONS
space_dimensions = "space.yml"
# RESULT DIR
result_path = "/home/iwe7/scfxn_optimize/results/"
pose_dir = "/home/iwe7/scfxn_optimize/results/poses/"
# OPTIMIZER args, higher values mean more exploration
xi = 0.01
kappa = 1.69
final_xi = 0.0001
final_kappa = 0.01
acq_func_kwargs = {"xi": xi, "kappa": kappa}
