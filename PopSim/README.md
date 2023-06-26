Multiscale Transmission model split across multiple scripts, each representing a different unit or scale-- Much of the code is inherited from https://github.com/InstituteforDiseaseModeling/community-structure-mediates-polio-transmission

Individual.py: Individual attributes

Household.py: Household attributes, tree structure, household splitting logic

Node.py: Nodes used in the Household tree structure

Bari.py: Bari level

Village.py Village level

Utility.py: collection of functions used throughout the simulation

Config.py: all the parameter values are placed here

Simulation.py: Village Simulation -- creates a set of serialized Village objects to run the epidemic on top of. Output pkl results are automatically stored to a folder labeled "sims". Will error if sims folder not found... will fix later



Requires:
ete3, scipy, numpy, collections, dill
