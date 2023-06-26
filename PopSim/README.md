Multiscale Transmission model split across multiple scripts, each representing a different unit or scale

Individual.py: Individual attributes

Household.py: Household attributes, tree structure, household splitting logic

Node.py: Nodes used in the Household tree structure

Bari.py: Bari level

Village.py Village level

Utility.py: collection of functions used throughout the simulation

Config.py: all the parameter values are placed here

Simulation.py: Village Simulation -- creates a set of serialized Village objects to run the epidemic on top of. Output pkl results are automatically stored to a folder labeled "sims". Will error if sims folder not found... will fix later

To run:
python Simulation.py

epidemic.py: Epidemic simulator: Takes 2 argument from the command line (for the moment). First argument: Requires a file with a list of serialized villages to run (sims/village_simulations_github_1.pkl for an example). 2nd argument: json file with the transmission parameters in it (epi_config.json for an example)

To run:
python epidemic.py {serialized_villages_file} {epi_config_json_file}

epi_config.json:
Overwrites the default config file with the desired parameters. overwritten parameters currently include: fecal_oral_dose (float, in the form 1e-#), hh_transmission (boolean flag), beta_hh (int), bari_transmission (boolean), beta_bari (int), village_transmission (boolean), beta_village (int), age_assortivity (boolean flag)

Requires:
ete3, scipy, numpy, collections, dill
