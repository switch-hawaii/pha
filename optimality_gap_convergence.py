"""User-defined convergence class to use with pysp's progressive hedging algorithm.

This depends on pysp's phboundextension. It can be loaded by specifying 
--user-defined-extension=pyomo.pysp.plugins.phboundextension --user-defined-extension=optimality_gap_convergence
on the runph command line.

By default, this will calculate inner and outer bounds every iteration, which is fairly expensive
(two solves with fixed investment variables). These calculations can be spread out by setting an 
environment variable PHBOUNDINTERVAL greater than 1.

This code is based on pyomo.pysp.convergence. PySP does not provide a user-definable convergence class, 
but we inject one by adding it to ph._convergers from the pre_ph_initialization method of a standard plugin.

"""

from pyomo.util.plugin import SingletonPlugin, implements
from pyomo.pysp import phextension
from pyomo.pysp.convergence import ConvergenceBase

convergence_threshold = 0.0001

class OptimalityGapConvergence(ConvergenceBase):

    def __init__(self, *args, **kwds):

        ConvergenceBase.__init__(self, *args, **kwds)
        self._name = "optimality-gap"

    def computeMetric(self, ph, scenario_tree, instances):
        # note: outer bound is the current solution, which may be inconsistent across scenarios.
        # for a minimization problem, this is the lower limit on the optimal objective
        # inner bound is the best available solution that is consisten across scenarios (xhat).
        # For a minimization problem, this is the upper limit on the optimal objective.
        
        # This calculates the size of the optimality gap relative to the best possible objective
        # value (the outer bound).
        if ph._reported_inner_bound is None or ph._reported_outer_bound is None:
            return None
        else:
            return abs(ph._reported_inner_bound - ph._reported_outer_bound) / ph._reported_outer_bound


class OptimalityPHExtension(SingletonPlugin):

    implements(phextension.IPHExtension) 

    def reset(self, ph):
        pass

    def pre_ph_initialization(self, ph):
        print "Adding OptimalityGapConvergence to list of convergers, with target={}.".format(convergence_threshold)
        # code here is similar to _converger-related code in pyomo.pysp.ph
        # note: this is applied in addition to the other convergers, so those 
        # should be set to strict values if we want this to be the main criterion
        # TODO: find a way to specify the threshold on the command-line
        # (a crude way would be to write several different versions of this file
        # with different hard-coded thresholds.)
        ph._convergers.append(OptimalityGapConvergence(convergence_threshold=convergence_threshold))

    def post_instance_creation(self, ph):
        pass

    def post_ph_initialization(self, ph):
        pass

    def post_iteration_0_solves(self, ph):
        pass

    def post_iteration_0(self, ph):
        pass

    def pre_iteration_k_solves(self, ph):
        pass

    def post_iteration_k_solves(self, ph):
        pass

    def post_iteration_k(self, ph):
        pass

    def post_ph_execution(self, ph):
        pass