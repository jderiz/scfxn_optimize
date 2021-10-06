"""
CONFIGURATION
"""
from relax import relax_with_config, initialize as relax_init
from design import design_with_config, initialize as design_init
from dummy import dummy_objective

### OBJECTIVE FUNCTION TO EVALUATE 
_dummy_objective = dummy_objective
_objective = relax_with_config
### SETUP THAT IS BEING DONE ONLY ONCE ON EACH CPU e.g initializing rosetta etc.
_init_method = relax_init 
### SPACE DIMENSIONS
space_dimensions = "space.yml"

#OPTIMIZER args, 
xi = 0.01
kappa = 1.69
final_xi = 0.0001 
final_kappa = 0.01
acq_func_kwargs = {"xi": xi, "kappa": kappa}
