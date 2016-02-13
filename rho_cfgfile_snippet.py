build_vars = [
    "BuildProj", "BuildBattery", 
    "BuildPumpedHydroMW", "BuildAnyPumpedHydro",
    "RFMSupplyTierActivate",
    "BuildElectrolyzerMW", "BuildLiquifierKgPerHour", "BuildLiquidHydrogenTankKg",
    "BuildFuelCellMW"
]

costs = []
baseval = value(m.Minimize_System_Cost)
# surprisingly slow, but it gets the job done
for var in build_vars:
    for v in getattr(m, var).values():
        # perturb the value of each variable to find its coefficient in the objective function
        v.value += 1; c = value(m.Minimize_System_Cost) - baseval; v.value -= 1
        costs.append((v.cname(), c))
