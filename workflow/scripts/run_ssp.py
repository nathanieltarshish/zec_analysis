
import utils
import pickle
import os

# Access parameters and wildcards from Snakemake
try:
    scenario = snakemake.params.scenario
except:
    scenario = snakemake.wildcards.scenario
    
final_year = snakemake.params.final_year
ZEC_year = snakemake.params.ZEC_year
output = snakemake.output[0]

# Generate the ensemble
f = utils.gen_fair_ensemble([scenario], final_year=final_year, ZEC_year=ZEC_year)
f.run()

# Ensure the output directory exists
output_pkl = output+".pkl"
os.makedirs(os.path.dirname(output), exist_ok=True)

# Save the results
with open(output_pkl, "wb") as file:
    pickle.dump(f, file)

output_nc = output+".nc"
f.to_netcdf(output_nc)
