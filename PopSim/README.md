Multiscale Transmission model split across multiple scripts, each representing a different unit or scale-- Much of the code is inherited from https://github.com/InstituteforDiseaseModeling/community-structure-mediates-polio-transmission

Individual.py: Individual attributes

Household.py: Household attributes, tree structure, household splitting logic

Node.py: Nodes used in the Household tree structure

Bari.py: Bari level

Village.py Village level

Utility.py: collection of functions used throughout the simulation

Config.py: all the parameter values are placed here

Simulation.py: Village Simulation -- creates a set of serialized Village objects to run the epidemic on top of. Output pkl results are automatically stored to a folder labeled "sims". Will error if sims folder not found... will fix later

To run:
Refer to the ipynb notebooks for model calibration and genetics.

Evo-Epi simulations
Make sure that a set of demographic files are created using Demographics.py. Use project_future_demography.py to project immunity by X years

PopSim/outbreak.py with a param file that describes the conditions to use. These can be created with the create_outbreak_configs.py and outbreak_configs.json is an example.
This produces three outputs, simulation results for the genetics, prevalence and immunity


Requires:
ete3, scipy, numpy, collections, dill
