
import utils
import pickle
import os

# Access parameters and wildcards from Snakemake

scenario = snakemake.wildcards.scenario
species = snakemake.wildcards.species
end = snakemake.wildcards.end

final_year = snakemake.params.final_year
output_pkl = snakemake.output.output_pkl
output_nc = snakemake.output.output_nc


if species == "all":
    f = utils.gen_fair_ensemble([scenario], final_year=final_year)
elif species == "co2":
    f = utils.gen_fair_co2_only([scenario], final_year=final_year)
elif species == "all-relax-land":
    f = utils.gen_fair_ensemble([scenario],
                                relax_land=True,
                                final_year=final_year,
                                relax_land_year=2024)    
    
if end == "net-zero-zec":
    net_zero_GHG_time = utils.get_net_zero_GHG_time(scenario)    
    utils.ZEC_emissions(f, net_zero_GHG_time)
elif end == "net-zero-decay":
    net_zero_GHG_time = utils.get_net_zero_GHG_time(scenario)    
    utils.decay_emissions(f, net_zero_GHG_time)    
elif end == "2100-zec":
    utils.ZEC_emissions(f, 2100)
if end == "2024-zec":
    utils.ZEC_emissions(f, 2024)    

f.run()

# Ensure the output directory exists
os.makedirs(os.path.dirname(output_pkl), exist_ok=True)

# Save the results
with open(output_pkl, "wb") as file:
    pickle.dump(f, file)

f.to_netcdf(output_nc)
