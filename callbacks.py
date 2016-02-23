#import debug

from pyomo.util.plugin import *
from pyomo.pysp import phextension
from pyomo.environ import value
from pyomo.pysp.phboundbase import ExtractInternalNodeSolutionsforInner
from pyomo.pysp.phutils import indexToString

import sys, os, datetime, inspect

cur_dir = os.path.abspath(os.path.dirname(__file__))
sys.path.append(os.path.join(cur_dir, "switch-hawaii-core"))
import util
from util import get

build_vars = [
    "BuildProj", "BuildBattery", 
    "BuildPumpedHydroMW", "BuildAnyPumpedHydro",
    "RFMSupplyTierActivate",
    "BuildElectrolyzerMW", "BuildLiquifierKgPerHour", "BuildLiquidHydrogenTankKg",
    "BuildFuelCellMW"
]
inputs_dir = "inputs"

# note: we make a best effort to get a unique, shared jobid for each job,
# but if the job isn't launched via mpi or slurm, different PHSolverServerExtensions may 
# get different job ids, since they are launched at different times.
jobid = os.environ.get('JOBID') # could be set by user
if jobid is None:   # not running under slurm
    jobid = os.environ.get('SLURM_JOBID')
if jobid is None:   # not running under slurm
    jobid = os.environ.get('OMPI_MCA_ess_base_jobid')
if jobid is None:
    jobid = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

def make_env_boolean(var):
    val = os.environ.get(var)
    if val is None:
        val = "False"
    if val.lower() in ["1", "true", "y", "yes", "on"]:
        return True
    elif val.lower() in ["0", "false", "n", "no", "off"]:
        return False
    else:
        raise ValueError('Unrecognized value for environment variable: {}="{}"'.format(var, val))
    
write_scenario_summary = make_env_boolean('PHA_WRITE_SCENARIO_SUMMARY')
write_iter0_scenario_builds = make_env_boolean('PHA_WRITE_ITER0_SCENARIO_BUILDS')
write_xhat_build = make_env_boolean('PHA_WRITE_XHAT_BUILD')

# callbacks for main runph process
# based on pyomo.pysp.plugins.testphextension and pyomo.pysp.plugins.examplephextension
class PHExtension(SingletonPlugin):

    implements(phextension.IPHExtension) 

    def reset(self, ph):
        pass

    def pre_ph_initialization(self, ph):
        pass

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
        #import pdb; pdb.Pdb(stdout=sys.__stdout__).set_trace() # have to grab original stdout, or this crashes

        if write_xhat_build:
            # note: when solved in serial mode, ph.get_scenario_tree().get_arbitrary_scenario()._instance
            # contains a solved instance that can (maybe) be used for reporting at this point,
            # as done in old_post_ph_execution below. However, when using pyro, that is not available.
            # Also, at this point with the parallel or serial solver, ph.get_scenario_tree()._solution
            # and ph._scenario_tree.findRootNode().get_variable_value(name, index) (which reads it) have no values.

            # So we take a different course, using the best available solution directly from the scenario tree.
            # The code below is modeled on pyomo.pysp.ph.ProgressiveHedging.pprint()
        
            # note: at this point, we might be able to use ph.get_scenario_tree().get_arbitrary_scenario()._x[root_node.name] 
            # to get a solution for one of the scenarios (the first), but ExtractInternalNodeSolutionsforInner() 
            # seems designed to give a definitive solution (i.e., the best admissable solution) (see 
            # compute_and_report_inner_bound_using_xhat which is calculated a few lines below the point 
            # where this is called in pyomo.pysp.ph)
        
            root_node = ph.get_scenario_tree().findRootNode()
            solution = ExtractInternalNodeSolutionsforInner(ph)[root_node.name] # {id: value}
            output_file = os.path.join("outputs", "build_{}_xhat.tsv".format(jobid))
            make_outputs_dir()
            write_build_vars_table(ph, solution, output_file)

# callbacks for phsolverservers
# based on pyomo.pysp.plugins.examplephextension
class PHSolverServerExtension(SingletonPlugin):

    implements (phextension.IPHSolverServerExtension)

    def pre_ph_initialization(self,ph):
        """Called before PH initialization."""

    def post_instance_creation(self,ph):
        """Called after the instances have been created."""

    def post_ph_initialization(self, ph):
        """Called after PH initialization."""

    def pre_iteration_0_solve(self, ph):
        """Called before the iteration 0 solve begins."""

    def post_iteration_0_solve(self, ph):
        """Called after the iteration 0 solve is finished."""
        # the caller knows which scenario this is, but they didn't tell us,
        # so we have to dig out the name (from phsolverserver solve() method)
        scenario_name = inspect.stack()[1][0].f_locals["object_name"]

        if write_iter0_scenario_builds:

            scenario = ph._scenario_tree._scenario_map[scenario_name] # based on phsolverserver.solve()
            root_node_name = ph._scenario_tree.findRootNode().name
            solution = scenario._x[root_node_name]  # based on tree_structure.scenario.update_solution_from_instance() 
                                                # and tree_structure.scenario.push_solution_to_instance
            # note: for some reason scenario._instance.solutions[0] is not an indexable object

            output_file = os.path.join("outputs", "build_{}_iter0_{}.tsv".format(jobid, scenario_name))
            make_outputs_dir()
            write_build_vars_table(ph, solution, output_file)


        if write_scenario_summary:
            write_summary_table(scenario_name, ph._instances[scenario_name])
        
    def pre_iteration_k_solve(self, ph):
        """Called before the iteration k solve begins."""

    def post_iteration_k_solve(self, ph):
        """Called after the iteration k solve is finished."""
        # note: we could write the results from all the scenarios into a single 
        # job file, but then we would risk write conflicts on Windows, which doesn't 
        # have easy file locking. It's also possible that the different solver
        # servers won't have the same jobid, which would make it difficult
        # to be sure of erasing the right file each iteration.

        # the caller knows which scenario this is, but they didn't tell us,
        # so we have to dig out the name (from phsolverserver solve() method)
        scenario_name = inspect.stack()[1][0].f_locals["object_name"]

        # print "post_iteration_k_solve for " + scenario_name

        if write_scenario_summary:
            write_summary_table(scenario_name, ph._instances[scenario_name])

def make_outputs_dir():
    if not os.path.isdir("outputs"):
        os.makedirs("outputs")

def write_build_vars_table(ph, solution, output_file):
    # root_node = ph.get_scenario_tree().findRootNode()
    # this works on runph and solverserver
    root_node = ph._scenario_tree.findRootNode()
    variable_names = sorted(root_node._variable_indices.keys())
    # alternatively, this could be (as in ph.pprint)
    # variable_names = sorted(ph.get_scenario_tree().stages[0]._variables.keys())

    # useful elements for getting variable values:
    # node._variable_indices = {varname: [(key1a, key1b), (key2a, key2b)]}
    # node._variable_ids = {varid: (varname, (key1a, key1b))}
    # node._name_index_to_id = {(varname, (key1a, key1b)): varid}

    variable_data = [
        (
            v + indexToString(k),   # var name with index
            solution[root_node._name_index_to_id[(v, k)]]   # current value
        )
            for v in variable_names 
                for k in sorted(root_node._variable_indices[v])
    ]
    
    print "writing {}...".format(output_file)
    # print "variables to write:"
    # print variable_data
    
    with open(output_file, 'w') as f:
        f.writelines(
            "\t".join(map(str, r)) + "\n"
                for r in [("variable", "value")] + variable_data
        )
                            
def write_summary_table(scenario_name, instance):
    m = instance

    # note: this gets overwritten after every solve. It would be more efficient to
    # write it only after the last solve, but there's no hook in phsolverserver for that
    # (it may not even be called at that point)
    make_outputs_dir()
    
    # the build variables may be pinned to match a particular file; 
    # if so, we use that as part of the filename
    build_file = os.environ.get('PHA_FIX_BUILD_VARS_FILE')
    if build_file is None:
        build_tag = ""
    else:
        # build_file typically looks like
        # "outputs/build_JJJJJJ_xhat.tsv" or "outputs/build_JJJJJJ_iter0_Scenario_nnnn.tsv"
        # we shorten this to tag like "bJJJJJJ_xhat" or "bJJJJJJ_iter0_Scenario_nnnn"
        # before prepending the current job id and appending the evaluated scenario ID.
        if os.path.dirname(build_file) == "outputs":
            # remove outputs dir if specified, otherwise keep dir name in there
            build_tag = os.path.basename(build_file)
        else:
            # remove path separators from the build file name
            build_tag = build_file.replace(os.sep, '_')
        if build_tag.startswith("build_"):
            build_tag = 'b' + build_tag[6:] # shorten "build_" to "b"
        if build_tag.endswith(".tsv"):
            build_tag = build_tag[:-4]
        build_tag = "_" + build_tag

        output_file = os.path.join("outputs", "summary_{}{}_{}.tsv".format(jobid, build_tag, scenario_name))
    
    print "writing {}...".format(output_file)
    
    period_duration = {
        pe:
        sum(m.tp_weight[tp] for tp in m.PERIOD_TPS[pe])
            for pe in m.PERIODS
    }
    
    values = []
    demand_components = [c for c in ('lz_demand_mw', 'DemandResponse', 'ChargeEVs') if hasattr(m, c)]
    
    #  total cost / kWh generated in each period 
    # (both discounted to today, so the discounting cancels out)
    values.extend([
        (
            "cost_per_kwh", pe,
            m.SystemCostPerPeriod[pe]
            / sum(
                m.bring_timepoint_costs_to_base_year[tp] * 1000.0 *
                sum(getattr(m, c)[lz, tp] for c in demand_components for lz in m.LOAD_ZONES)
                    for tp in m.PERIOD_TPS[pe]
            )
        )
            for pe in m.PERIODS
    ])
    # Renewable energy share
    if hasattr(m, 'RPSEligiblePower'):
        # total renewable share over all periods
        values.extend([
            (
                "renewable_share", pe, m.RPSEligiblePower[pe]/m.RPSTotalPower[pe]
            )
                for pe in m.PERIODS
        ])
    # average production from each fuel during each period
    values.extend([
        (
            f, pe,
            sum(
                get(m.DispatchProjByFuel, (pr, tp, f), 0.0) * m.tp_weight[tp]
                for pr in m.PROJECTS_BY_FUEL[f] for tp in m.PERIOD_TPS[pe]
            ) / period_duration[pe]
        )
            for f in m.FUELS for pe in m.PERIODS
    ])    
    # total production from each non-fuel source
    values.extend([
        (
            s, pe,
            sum(
                get(m.DispatchProj, (pr, tp), 0.0) * m.tp_weight[tp]
                for pr in m.PROJECTS_BY_NON_FUEL_ENERGY_SOURCE[s] for tp in m.PERIOD_TPS[pe]
            ) / period_duration[pe]
        )
            for s in m.NON_FUEL_ENERGY_SOURCES for pe in m.PERIODS
    ])    
    # curtailments
    values.extend([
        (
            "curtail_"+s, pe,
            sum(
                (get(m.DispatchUpperLimit, (pr, tp), 0.0) - get(m.DispatchProj, (pr, tp), 0.0)) 
                * m.tp_weight[tp]
                for pr in m.PROJECTS_BY_NON_FUEL_ENERGY_SOURCE[s] for tp in m.PERIOD_TPS[pe]
            ) / period_duration[pe]
        )
            for s in m.NON_FUEL_ENERGY_SOURCES for pe in m.PERIODS
    ])    
    # all LZ_Energy_Components
    values.extend([
        (
            component, pe,
            sum(
                getattr(m, component)[lz, tp] * m.tp_weight[tp] 
                    for lz in m.LOAD_ZONES for tp in m.PERIOD_TPS[pe]
            ) / period_duration[pe]
        )
            for component in m.LZ_Energy_Components_Produce for pe in m.PERIODS
    ])    
    values.extend([
        (
            component, pe,
            sum(
                getattr(m, component)[lz, tp] * m.tp_weight[tp] 
                    for lz in m.LOAD_ZONES for tp in m.PERIOD_TPS[pe]
            ) / period_duration[pe]
        )
            for component in m.LZ_Energy_Components_Consume for pe in m.PERIODS
    ])    

    with open(output_file, 'w') as f:
        f.writelines(
            "\t".join((key, str(per), str(value(val)))) + "\n"
                for (key, per, val) in values
        )



rho_cost_multiplier = 1.0   # default value

rhos = None

# based on pyomo_examples_11103/pysp/sizes/config/rhosetter.py and pyomo/pysp/ph.py
def ph_rhosetter_callback(ph, scenario_tree, scenario):
    global rhos
   
    m = scenario._instance
    
    if rhos is None:    # read rho values from disk if not cached already
        # read previously stored rho values (much faster than calculating them)
        with open(os.path.join(inputs_dir, "rhos.tsv"), "r") as f:
            rows = [r[:-1].split('\t') for r in f]  # split at tabs; omit newlines
            rhos = {r[0]: float(r[1]) for r in rows}

    for var_name in build_vars:
        vars = getattr(m, var_name).values()
        for v in vars:
            ph.setRhoOneScenario(
                scenario_tree.findRootNode(),
                scenario,
                m._ScenarioTreeSymbolMap.getSymbol(v),
                max(1e-6, rho_cost_multiplier * rhos[v.cname()]))   # never set rho to zero, even if var isn't in obj. fn.
