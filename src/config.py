"""
CONFIGURATION
Import your objective function. If you have an initializer function to run on 
worker instantiation also import it here and assign it to the correspinging 
variables.
Set custom for xi and kappa if you whish else leave them at the default.
Specify a path where your results should get saved to 
"""
import logging

from colorlog import ColoredFormatter

from allosteric import initialize as relax_init
from allosteric import relax_with_config
from design import design_with_config
from design import initialize as design_init
from dummy import dummy_objective

# LOG LEVEL
level = logging.DEBUG
# formatter = ColoredFormatter(
#    '%(log_color)s%(levelname)s:-PID%(process)d[%(threadName)s]::%(name)s:%(funcName)s:%(white)s%(message)s',
#    datefmt=None,
#    reset=True,
#    log_colors={
#        'DEBUG':    'cyan',
#        'INFO':     'green',
#        'WARNING':  'yellow',
#        'ERROR':    'red',
#        'CRITICAL': 'red,bg_white',
#    },
#    secondary_log_colors={},
#    style='%'
# )
#color_handler = logging.StreamHandler()
# color_handler.setFormatter(formatter)
# logging.getLogger().addHandler(color_handler)
logging.getLogger().setLevel(level)
# OBJECTIVE FUNCTION TO EVALUATE
_dummy_objective = dummy_objective
_objective = relax_with_config
# SETUP THAT IS BEING DONE ONLY ONCE ON EACH CPU e.g initializing rosetta etc.
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
