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

emissions = emissions.sel(
    scenario=["ssp119", "ssp126", "ssp245", "ssp534-over"]
)
emissions = emissions.assign_coords(
    scenario=["ssp119", "ssp126", "historical", "ssp534-over"]
)

emissions.loc[
    dict(timepoints=slice(2025, None), scenario="historical")
] = emissions.loc[
    dict(timepoints=2024.5, scenario="historical")
]
emissions.to_netcdf(utils.get_path("results/emissions.nc"))
