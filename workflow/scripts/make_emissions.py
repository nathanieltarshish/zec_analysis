#!/usr/bin/env python
import utils
import xarray as xr

prefix = utils.get_path("datasets/fair_calibrate")
fair_v = "2.1.3"
cal_v = "1.4"
constraint_set = "all-2022"

emissions_file = (
    prefix + f"/output/fair-{fair_v}/v{cal_v}/{constraint_set}/emissions/ssps_harmonized_1750-2499.nc"
)
emissions = xr.load_dataarray(emissions_file)

historical = emissions.sel(
    scenario=["ssp245"]
)
historical = historical.assign_coords(
    scenario=["historical"]
)

historical.loc[
    dict(timepoints=slice(2025, None), scenario="historical")
] = historical.loc[
    dict(timepoints=2024.5, scenario="historical")
]

emissions_final = xr.concat([emissions, historical], dim="scenario")

emissions_final.to_netcdf(utils.get_path("results/emissions.nc"))
