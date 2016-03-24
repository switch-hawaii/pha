from pyomo.environ import *

def define_components(m):
    ###################
    # resource rules to match HECO's PSIP (8/26/15, http://files.hawaii.gov/puc/3_Dkt%202011-0206%202014-08-26%20HECO%20PSIP%20Report.pdf)
    ##################
    
    # run AES on 50% coal, 50% pellet biomass
    m.AES_Half_Biomass = Constraint(m.TIMEPOINTS, rule=lambda m, tp:
        m.ProjFuelUseRate["Oahu_AES_GEN1", tp, "Pellet-Biomass"] 
        == m.ProjFuelUseRate["Oahu_AES_GEN1", tp, "Coal"]
    )
    
    # get 121.0 GWh/year from biodiesel
    m.Biodiesel_Target = Constraint(m.PERIODS, rule=lambda m, per:
        sum(
            m.DispatchProjByFuel[p, tp, "Biodiesel"] * m.tp_weight[tp]
                for p in m.PROJECTS_BY_FUEL["Biodiesel"]
                    for tp in m.PERIOD_TPS[per]
                        if (p, tp) in m.PROJ_DISPATCH_POINTS
        )/m.period_length_years[per] == 121000.0
    )
    
    # targets for individual generation technologies, extrapolated into the future
    technology_target_gwh_2021 = {
        "DistPV_peak": 820.1,
        "CentralTrackingPV": 501.8,
        "Wind": 319.6,
    }
    technology_target_gwh_2030 = {
        "DistPV_peak": 1031.6,
        "CentralTrackingPV": 501.8,
        "Wind": 460.6,
    }
    m.TARGETED_TECHNOLOGIES = Set(initialize=lambda m: technology_target_gwh_2021.keys())
    m.technology_target_gwh = Param(m.TARGETED_TECHNOLOGIES, m.PERIODS, rule = lambda m, tech, per:
        technology_target_gwh_2021[tech] 
        + (m.per-2021.0) * (technology_target_gwh_2030[tech] - technology_target_gwh_2021[tech]) 
            / (2030.0-2021.0)
    )
    
    m.Technology_Targets = Constraint(m.TARGETED_TECHNOLOGIES, m.PERIODS, rule=lambda m, tech, per:
        sum(
            m.DispatchProj[p, tp] * m.tp_weight[tp] / m.period_length_years[per]
                for p in m.PROJECTS if m.proj_gen_tech[p] == tech
                    for tp in m.PERIOD_TPS[per] if (p, tp) in m.PROJ_DISPATCH_POINTS
        ) == m.technology_target_gwh[tech, per]
    )

