#!/usr/bin/env python
""" Utility functions used throughout the analysis."""

import re
import numpy as np
import xarray as xr
import pandas as pd
from fair import FAIR
from fair.interface import fill, initialise
from fair.io import read_properties
import os

def get_path(path):
    absolute_path = "/home/tarshish/Documents/research/zec/"
    return absolute_path + path 

def decay_emissions(f, year, tau=100):
    net_zero_emissions = f.emissions.sel(timepoints=year)
    future = f.emissions.where(f.emissions.timepoints > year, drop=True)                                                                                           
    t = future.timepoints - year
    baseline = f.species_configs.baseline_emissions
    decaying_emissions = net_zero_emissions*np.exp(-t/tau) + (1 - np.exp(-t/tau))*baseline    
    decaying_emissions = decaying_emissions.transpose("timepoints", "scenario", "config", "specie")
    past = f.emissions.where(f.emissions.timepoints <= year, drop=True)                                                                                            
    f.emissions = xr.concat([past, decaying_emissions], dim="timepoints")
    
def ZEC_emissions(f, year):
    future = f.emissions.where(f.emissions.timepoints > year, drop=True)                                                                                           
    future = future.fillna(0.0) * 0 + f.species_configs.baseline_emissions                                                                                             
    past = f.emissions.where(f.emissions.timepoints <= year, drop=True)                                                                                            
    f.emissions = xr.concat([past, future], dim="timepoints")

def gen_fair_ensemble(
    scenarios,
    final_year=2040,
    ensemble_members=841,
    relax_land=False,
    relax_land_year=None,
):

    prefix = get_path("datasets/fair_calibrate")
    fair_v = "2.1.3"
    cal_v = "1.4"
    constraint_set = "all-2022"

    calibration_dir=prefix+f"/output/fair-{fair_v}/v{cal_v}/{constraint_set}/calibrations/"
    posterior_dir=prefix+f"/output/fair-{fair_v}/v{cal_v}/{constraint_set}/posteriors/"

    methane_config = calibration_dir+"CH4_lifetime.csv"
    calibrated_params = posterior_dir+"calibrated_constrained_parameters.csv"
    landuse_factors = calibration_dir+"landuse_scale_factor.csv"
    lapsi_factors = calibration_dir+"lapsi_scale_factor.csv"

    # Loading the configuration dataframes
    df_methane = pd.read_csv(methane_config, index_col=0)
    df_configs = pd.read_csv(calibrated_params, index_col=0)
    df_landuse = pd.read_csv(landuse_factors, index_col=0)
    df_lapsi = pd.read_csv(lapsi_factors, index_col=0)

    da_emissions = xr.load_dataarray(get_path("results/emissions.nc"))

    # instantiate model with SSP defaults
    df_configs_subset = df_configs[:ensemble_members]
    valid_all = df_configs_subset.index
    f = FAIR(
        ch4_method="Thornhill2021",
        relax_land=relax_land,
        relax_land_year=relax_land_year,
    )
    f.define_time(1750, final_year, 1)
    f.define_scenarios(scenarios)
    f.define_configs(valid_all)
    species, properties = read_properties()

    # removed natural and species not accounted for in FaIR version and natural forcings 
    species.remove("Halon-1202")
    species.remove("NOx aviation")
    species.remove("Contrails")
    species.remove("Volcanic")
    species.remove("Solar")
    del properties["Volcanic"]
    del properties["Solar"]
    f.define_species(species, properties)
    f.allocate()
    f.fill_species_configs()
    
    # climate response
    fill(
        f.climate_configs["ocean_heat_capacity"],
        df_configs_subset.loc[:, "clim_c1":"clim_c3"].values,
    )
    fill(
        f.climate_configs["ocean_heat_transfer"],
        df_configs_subset.loc[:, "clim_kappa1":"clim_kappa3"].values,
    )  # not massively robust, since relies on kappa1, kappa2, kappa3 being in adjacent cols
    fill(
        f.climate_configs["deep_ocean_efficacy"],
        df_configs_subset["clim_epsilon"].values.squeeze(),
    )
    fill(
        f.climate_configs["gamma_autocorrelation"],
        df_configs_subset["clim_gamma"].values.squeeze(),
    )
    fill(
        f.climate_configs["sigma_eta"],
        df_configs_subset["clim_sigma_eta"].values.squeeze(),
    )
    fill(
        f.climate_configs["sigma_xi"],
        df_configs_subset["clim_sigma_xi"].values.squeeze(),
    )
    fill(f.climate_configs["seed"], df_configs_subset["seed"])
    fill(f.climate_configs["stochastic_run"], False)
    fill(f.climate_configs["use_seed"], False) 
    fill(f.climate_configs["forcing_4co2"], df_configs_subset["clim_F_4xCO2"])

    # carbon cycle
    fill(
        f.species_configs["iirf_0"],
        df_configs_subset["cc_r0"].values.squeeze(),
        specie="CO2",
    )
    fill(
        f.species_configs["iirf_airborne"],
        df_configs_subset["cc_rA"].values.squeeze(),
        specie="CO2",
    )
    fill(
        f.species_configs["iirf_uptake"],
        df_configs_subset["cc_rU"].values.squeeze(),
        specie="CO2",
    )
    fill(
        f.species_configs["iirf_temperature"],
        df_configs_subset["cc_rT"].values.squeeze(),
        specie="CO2",
    )

    # aerosol indirect
    fill(f.species_configs["aci_scale"], df_configs_subset["aci_beta"].values.squeeze())
    fill(
        f.species_configs["aci_shape"],
        df_configs_subset["aci_shape_so2"].values.squeeze(),
        specie="Sulfur",
    )
    fill(
        f.species_configs["aci_shape"],
        df_configs_subset["aci_shape_bc"].values.squeeze(),
        specie="BC",
    )
    fill(
        f.species_configs["aci_shape"],
        df_configs_subset["aci_shape_oc"].values.squeeze(),
        specie="OC",
    )

    # methane lifetime baseline and sensitivity
    fill(
        f.species_configs["unperturbed_lifetime"],
        df_methane.loc["historical_best", "base"],
        specie="CH4",
    )
    fill(
        f.species_configs["ch4_lifetime_chemical_sensitivity"],
        df_methane.loc["historical_best", "CH4"],
        specie="CH4",
    )
    fill(
        f.species_configs["ch4_lifetime_chemical_sensitivity"],
        df_methane.loc["historical_best", "N2O"],
        specie="N2O",
    )
    fill(
        f.species_configs["ch4_lifetime_chemical_sensitivity"],
        df_methane.loc["historical_best", "VOC"],
        specie="VOC",
    )
    fill(
        f.species_configs["ch4_lifetime_chemical_sensitivity"],
        df_methane.loc["historical_best", "NOx"],
        specie="NOx",
    )
    fill(
        f.species_configs["ch4_lifetime_chemical_sensitivity"],
        df_methane.loc["historical_best", "HC"],
        specie="Equivalent effective stratospheric chlorine",
    )
    fill(
        f.species_configs["lifetime_temperature_sensitivity"],
        df_methane.loc["historical_best", "temp"],
    )

    # correct land use and LAPSI scale factor terms
    fill(
        f.species_configs["land_use_cumulative_emissions_to_forcing"],
        df_landuse.loc["historical_best", "CO2_AFOLU"],
        specie="CO2 AFOLU",
    )
    fill(
        f.species_configs["lapsi_radiative_efficiency"],
        df_lapsi.loc["historical_best", "BC"],
        specie="BC",
    )

    # baseline emissions at 1750 in emissions dataset
    fill(f.species_configs["baseline_emissions"], 38.246272, specie="CH4")
    fill(f.species_configs["baseline_emissions"], 1.000925, specie="N2O")
    fill(f.species_configs["baseline_emissions"], 19.41683292, specie="NOx")
    fill(f.species_configs["baseline_emissions"], 2.293964929, specie="Sulfur")
    fill(f.species_configs["baseline_emissions"], 348.4549732, specie="CO")
    fill(f.species_configs["baseline_emissions"], 60.62284009, specie="VOC")
    fill(f.species_configs["baseline_emissions"], 2.096765609, specie="BC")
    fill(f.species_configs["baseline_emissions"], 15.44571911, specie="OC")
    fill(f.species_configs["baseline_emissions"], 6.656462698, specie="NH3")
    fill(f.species_configs["baseline_emissions"], 0.02125667, specie="CCl4")
    fill(f.species_configs["baseline_emissions"], 204.42373382567882, specie="CHCl3")
    fill(f.species_configs["baseline_emissions"], 207.4912295882706, specie="CH2Cl2")
    fill(f.species_configs["baseline_emissions"], 4407.5984529412835, specie="CH3Cl")
    fill(f.species_configs["baseline_emissions"], 106.68917253105691, specie="CH3Br")
    fill(
        f.species_configs["baseline_emissions"],
        0.008275907538039383,
        specie="Halon-1211",
    )
    fill(
        f.species_configs["baseline_emissions"], 1.0386253894514322e-05, specie="SO2F2"
    )
    fill(f.species_configs["baseline_emissions"], 0, specie="CF4")

    # aerosol direct
    for specie in [
        "BC",
        "CH4",
        "N2O",
        "NH3",
        "NOx",
        "OC",
        "Sulfur",
        "VOC",
        "Equivalent effective stratospheric chlorine",
    ]:
        fill(
            f.species_configs["erfari_radiative_efficiency"],
            df_configs_subset[f"ari_{specie}"],
            specie=specie,
        )

    # forcing scaling
    for specie in [
        "CO2",
        "CH4",
        "N2O",
        "Stratospheric water vapour",
        "Light absorbing particles on snow and ice",
        "Land use",
    ]:
        fill(
            f.species_configs["forcing_scale"],
            df_configs_subset[f"fscale_{specie}"].values.squeeze(),
            specie=specie,
        )

    for specie in [
        "CFC-11",
        "CFC-12",
        "CFC-113",
        "CFC-114",
        "CFC-115",
        "HCFC-22",
        "HCFC-141b",
        "HCFC-142b",
        "CCl4",
        "CHCl3",
        "CH2Cl2",
        "CH3Cl",
        "CH3CCl3",
        "CH3Br",
        "Halon-1211",
        "Halon-1301",
        "Halon-2402",
        "CF4",
        "C2F6",
        "C3F8",
        "c-C4F8",
        "C4F10",
        "C5F12",
        "C6F14",
        "C7F16",
        "C8F18",
        "NF3",
        "SF6",
        "SO2F2",
        "HFC-125",
        "HFC-134a",
        "HFC-143a",
        "HFC-152a",
        "HFC-227ea",
        "HFC-23",
        "HFC-236fa",
        "HFC-245fa",
        "HFC-32",
        "HFC-365mfc",
        "HFC-4310mee",
    ]:
        fill(
            f.species_configs["forcing_scale"],
            df_configs_subset["fscale_minorGHG"].values.squeeze(),
            specie=specie,
        )

    # ozone
    for specie in [
        "CH4",
        "N2O",
        "Equivalent effective stratospheric chlorine",
        "CO",
        "VOC",
        "NOx",
    ]:
        fill(
            f.species_configs["ozone_radiative_efficiency"],
            df_configs_subset[f"o3_{specie}"],
            specie=specie,
        )

    # initial condition of CO2 concentration(but not baseline for forcing calculations)
    fill(
        f.species_configs["baseline_concentration"],
        df_configs_subset["cc_co2_concentration_1750"].values.squeeze(),
        specie="CO2",
    )

    # initial conditions
    initialise(f.concentration, f.species_configs["baseline_concentration"])
    initialise(f.forcing, 0)
    initialise(f.temperature, 0)
    initialise(f.cumulative_emissions, 0)
    initialise(f.airborne_emissions, 0)

    #fill emissions
    # Emissions loading and initialization
    anthro_emis = da_emissions.where(
        ~da_emissions.specie.isin(["Solar", "Volcanic"]), drop=True
    )
    anthro_emis = anthro_emis.sel(scenario=scenarios)
    anthro_emis = anthro_emis.interp(timepoints=f.timepoints, method="nearest")
    anthro_emis = anthro_emis.loc[dict(config="unspecified")]
    f.emissions = (
        anthro_emis.expand_dims(dim=["config"], axis=(2))
        .drop("config")
        * np.ones((1, 1, ensemble_members, 1))
    )

    return f

def plot_median(
    x,
    y,
    axis,
    color="black",
    lw=1,
    alpha=1,
    zorder=5,
    spread=None,
    spread_color=None,
    spread_lw=None,
    spread_alpha=0.4,
    spread_ec=None,
    spread_zorder=4,
    label=None,
):
    """Plots the median and desired interquantile range (spread for a data timeseries"""

    axis.plot(
        x,
        np.median(y, axis=1),
        color=color,
        alpha=alpha,
        zorder=zorder,
        label=label,
        linewidth=lw,
    )

    if spread is not None:
        if spread_color is None:
            spread_color = color
        axis.fill_between(
            x,
            np.percentile(y, spread[0], axis=1),
            np.percentile(y, spread[1], axis=1),
            color=spread_color,
            alpha=spread_alpha,
            zorder=spread_zorder,
            linewidth=spread_lw,
            ec=spread_ec,
        )


def replace_chemical_notation(text):
    """Replace patterns like CO2, CH4, N2O, NH3 with their LaTeX subscripts."""
    text = re.sub(r"([A-Z][a-z]*)(\d+)", r"\1$_\2$", text)
    return text


def set_plot_configs(plt, fsize=10):
    """Set the matplotlib parameters to the dersired Tex and font configs."""
    plt.rcParams.update({
        'text.usetex': True,               # Enable LaTeX for all text rendering
        'font.family': 'sans-serif',       # Set font family to sans-serif
        'font.sans-serif': ['Helvetica'],  # Use Helvetica as the sans-serif font
        'text.latex.preamble': r'\usepackage{amsmath} \usepackage{helvet} \usepackage{sfmath} \renewcommand{\familydefault}{\sfdefault}',  
        'font.size': fsize                    # Set the default font size
    })


def sum_forcing_set(forcing, species, name):
    """Sum over a set of species and add it to the forcing dataset."""
    new_forcing = forcing.sel(specie=species).sum("specie")
    new_forcing = new_forcing.expand_dims(specie=[name], axis=-1)
    forcing = xr.merge([forcing, new_forcing])
    return forcing


def aggregate_forcings(forcing):
    """Group forcings into categories and them to the forcing dataset."""

    forcing = sum_forcing_set(forcing, forcing.specie, "Anthropogenic")

    long_lived = [
        "CO2",
        "N2O",
        "CFC-12",
        "CFC-114",
        "CFC-115",
        "CF4",
        "C2F6",
        "C3F8",
        "c-C4F8",
        "C4F10",
        "C5F12",
        "C6F14",
        "C7F16",
        "C8F18",
        "NF3",
        "SF6",
        "HFC-23",
        "HFC-236fa",
    ]

    # Halon-1202 is not included in Fair calibrate v1.4.1
    short_lived = [
        "CH4",
        "CFC-11",
        "CFC-113",
        "HCFC-22",
        "HCFC-141b",
        "HCFC-142b",
        "CCl4",
        "CHCl3",
        "CH2Cl2",
        "CH3Cl",
        "CH3CCl3",
        "CH3Br",
        "Halon-1211",
        "Halon-1301",
        "Halon-2402",
        "SO2F2",
        "HFC-125",
        "HFC-134a",
        "HFC-143a",
        "HFC-152a",
        "HFC-227ea",
        "HFC-245fa",
        "HFC-32",
        "HFC-365mfc",
        "HFC-4310mee",
        "Aerosol-radiation interactions",
        "Aerosol-cloud interactions",
        "Ozone",
        "Light absorbing particles on snow and ice",
        "Stratospheric water vapour",
        "Land use",
        "Equivalent effective stratospheric chlorine",
    ]

    halogenated = [
        "CFC-11",
        "CFC-12",
        "CFC-113",
        "CFC-114",
        "CFC-115",
        "HCFC-22",
        "HCFC-141b",
        "HCFC-142b",
        "CCl4",
        "CHCl3",
        "CH2Cl2",
        "CH3Cl",
        "CH3CCl3",
        "CH3Br",
        "Halon-1211",
        "Halon-1301",
        "Halon-2402",
        "CF4",
        "C2F6",
        "C3F8",
        "c-C4F8",
        "C4F10",
        "C5F12",
        "C6F14",
        "C7F16",
        "C8F18",
        "NF3",
        "SF6",
        "SO2F2",
        "HFC-125",
        "HFC-134a",
        "HFC-143a",
        "HFC-152a",
        "HFC-227ea",
        "HFC-23",
        "HFC-236fa",
        "HFC-245fa",
        "HFC-32",
        "HFC-365mfc",
        "HFC-4310mee",
    ]

    # Halogenated gases with a lifetime of less than 100 years
    short_halogens = [
        "CFC-113",
        "HCFC-22",
        "HCFC-141b",
        "HCFC-142b",
        "CCl4",
        "CHCl3",
        "CH2Cl2",
        "CH3Cl",
        "CH3CCl3",
        "CH3Br",
        "Halon-1211",
        "Halon-1301",
        "Halon-2402",
        "SO2F2",
        "HFC-125",
        "HFC-134a",
        "HFC-143a",
        "HFC-152a",
        "HFC-227ea",
        "HFC-245fa",
        "HFC-32",
        "HFC-365mfc",
        "HFC-4310mee",
    ]

    # Halogenated gases with a lifetime of 100 years or greater
    long_halogens = [
        "CFC-12",
        "CFC-114",
        "CFC-115",
        "CF4",
        "C2F6",
        "C3F8",
        "c-C4F8",
        "C4F10",
        "C5F12",
        "C6F14",
        "C7F16",
        "C8F18",
        "NF3",
        "SF6",
        "HFC-23",
        "HFC-236fa",
    ]

    miscellaneous = [
        "Light absorbing particles on snow and ice",
        "Stratospheric water vapour",
        "Land use",
        "Equivalent effective stratospheric chlorine",
    ]

    aerosols = ["Aerosol-radiation interactions", "Aerosol-cloud interactions"]

    forcing = sum_forcing_set(forcing, long_lived, "Long-lived")
    forcing = sum_forcing_set(forcing, short_lived, "Short-lived")
    forcing = sum_forcing_set(forcing, halogenated, "Halogens")
    forcing = sum_forcing_set(forcing, miscellaneous, "Miscellaneous")
    forcing = sum_forcing_set(forcing, aerosols, "Aerosols")
    forcing = sum_forcing_set(forcing, short_halogens, "Short halogens")
    forcing = sum_forcing_set(forcing, long_halogens, "Long halogens")

    forcing["specie"] = forcing["specie"].where(
        forcing["specie"] != "Land use", "Land-use albedo"
    )
    forcing["specie"] = forcing["specie"].where(
        forcing["specie"] != "Light absorbing particles on snow and ice",
        "Particles on \n snow + ice",
    )
    forcing["specie"] = forcing["specie"].where(
        forcing["specie"] != "Aerosol-cloud interactions", "Aerosol-cloud"
    )
    forcing["specie"] = forcing["specie"].where(
        forcing["specie"] != "Aerosol-radiation interactions", "Aerosol-radiation"
    )
    forcing["specie"] = forcing["specie"].where(
        forcing["specie"] != "Stratospheric water vapour", "Stratospheric water vapor"
    )

    return forcing.forcing

def get_color(scenario):
    ar6_colors = {
        "ssp119": "#00a9cf",
        "ssp126": "#003466",
        "ssp245": "#f69320",
        "ssp370": "#df0000",
        "ssp434": "#2274ae",
        "ssp460": "#b0724e",
        "ssp534-over": "#92397a",
        "ssp585": "#980002",
        "historical": "black",        
    }

    return ar6_colors[scenario]

def get_title(scenario):
    
   fancy_titles = {
       "ssp119": "SSP1-1.9",
       "ssp126": "SSP1-2.6",
       "ssp245": "SSP2-4.5",
       "ssp370": "SSP3-7.0",
       "ssp434": "SSP4-3.4",
       "ssp460": "SSP4-6.0",
       "ssp534-over": "SSP5-3.4-overshoot",
       "ssp585": "SSP5-8.5",
       "historical": "black",       
   }

   return fancy_titles[scenario]

def get_net_zero_GHG_time(scenario):
    
    net_zero_times = {'historical': 2025.5, 'ssp119': 2081.5, 'ssp126': 2088.5, 'ssp534-over': 2073.5}

    return net_zero_times[scenario]
