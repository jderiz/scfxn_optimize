"""
CONFIGURATION
"""
import logging

from colorlog import ColoredFormatter

from design import design_with_config
from design import initialize as design_init
from dummy import dummy_objective
from relax import initialize as relax_init
from relax import relax_with_config

formatter = ColoredFormatter(
    '%(log_color)s%(levelname)s:-PID%(process)d[%(threadName)s]::%(name)s:%(funcName)s:%(white)s%(message)s',
    datefmt=None,
    reset=True,
    log_colors={
        'DEBUG':    'cyan',
        'INFO':     'green',
        'WARNING':  'yellow',
        'ERROR':    'red',
        'CRITICAL': 'red,bg_white',
    },
    secondary_log_colors={},
    style='%'
)
color_handler = logging.StreamHandler()
color_handler.setFormatter(formatter)
logging.getLogger().addHandler(color_handler)
# logging.basicConfig(
#     format='%(levelname)s:-PID%(process)d[%(threadName)s]::%(name)s:%(funcName)s:%(message)s', level=logging.DEBUG)
# coloredlogs.install(fmt='%(levelname)s:-PID%(process)d[%(threadName)s]::%(name)s:%(funcName)s:%(message)s', level="DEBUG")
RAY_DASHBOARD_IP = '172.22.180.238:6379'
# OBJECTIVE FUNCTION TO EVALUATE
_dummy_objective = dummy_objective
_objective = relax_with_config
# SETUP THAT IS BEING DONE ONLY ONCE ON EACH CPU e.g initializing rosetta etc.
_init_method = relax_init
# SPACE DIMENSIONS
space_dimensions = "space.yml"
# Directory where results should get saved
result_path="../results"
# OPTIMIZER args,
xi = 0.01
kappa = 1.69
final_xi = 0.0001
final_kappa = 0.01
acq_func_kwargs = {"xi": xi, "kappa": kappa}
