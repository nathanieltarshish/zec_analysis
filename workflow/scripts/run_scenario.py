
import utils
import pickle
import os

# Access parameters and wildcards from Snakemake
try:
    scenario = snakemake.params.scenario
except:
    scenario = snakemake.wildcards.scenario
    
final_year = snakemake.params.final_year
output_pkl = snakemake.output.output_pkl
output_nc = snakemake.output.output_nc
end = snakemake.wildcards.end

net_zero_GHG_time = utils.get_net_zero_GHG_time(scenario)
    
# Generate the ensemble
f = utils.gen_fair_ensemble([scenario],
                            final_year=final_year)
if end == "zec":
    utils.ZEC_emissions(f, net_zero_GHG_time)
elif end == "decay":
    utils.decay_emissions(f, net_zero_GHG_time)    

f.run()

# Ensure the output directory exists
os.makedirs(os.path.dirname(output_pkl), exist_ok=True)

# Save the results
with open(output_pkl, "wb") as file:
    pickle.dump(f, file)

f.to_netcdf(output_nc)
