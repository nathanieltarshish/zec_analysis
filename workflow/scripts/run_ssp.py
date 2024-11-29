
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
output_pkl = snakemake.output.output_pkl
output_nc = snakemake.output.output_nc
end = snakemake.wildcards.end

if end == "netzero":
    ZEC_year = None
    zero_CO2_only=True
elif end == "ZEC":
    zero_CO2_only=False
    
# Generate the ensemble
f = utils.gen_fair_ensemble([scenario],
                            final_year=final_year,
                            ZEC_year=ZEC_year,
                            zero_CO2_only=zero_CO2_only)
f.run()

# Ensure the output directory exists
os.makedirs(os.path.dirname(output_pkl), exist_ok=True)

# Save the results
with open(output_pkl, "wb") as file:
    pickle.dump(f, file)

f.to_netcdf(output_nc)
