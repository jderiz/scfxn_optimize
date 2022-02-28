"""
CONFIGURATION
Import your objective function. If you have an initializer function to run on 
worker instantiation also import it here and assign it to the corresponding
variables.
Set custom for xi and kappa if you whish else leave them at the default.
Specify a path where your results should get saved to.

"""
import logging

from colorlog import ColoredFormatter

from allosteric import initialize as relax_init
from allosteric import relax_with_config
from design import design_with_config
from design import initialize as design_init

# LOG LEVEL
level = logging.DEBUG
handler = logging.StreamHandler()
logging.getLogger().addHandler(handler)
logging.getLogger().setLevel(level)
# OBJECTIVE FUNCTION TO EVALUATE
_objective = relax_with_config
# SETUP THAT IS BEING DONE ONLY ONCE ON EACH CPU e.g loading large files
_init_method = relax_init
# SPACE DIMENSIONS
space_dimensions = "space.yml"
# RESULT DIR
result_path = "../results/"
# OPTIMIZER args, higher values mean more exploration
xi = 0.01
kappa = 1.69
final_xi = 0.0001
final_kappa = 0.01
acq_func_kwargs = {"xi": xi, "kappa": kappa}
