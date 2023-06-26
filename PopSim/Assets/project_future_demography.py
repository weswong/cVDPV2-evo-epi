import os
import sys
from collections import defaultdict
from Demographics import Demographics, write_dill
import dill
import config
import json
import copy
import numpy as np
from Infection import ImmunoInfection
import copy
from Individual import Individual
from Household import Household


def search_files(directory='.', extension='.pkl'):
    found_files = []
    extension = extension.lower()
    for dirpath, dirnames, files in os.walk(directory):
        for name in files:
            if extension and name.lower().endswith(extension):
                found_files.append(os.path.join(dirpath, name))

    return found_files


dir_path = os.path.dirname(os.path.realpath(__file__))
output_dir = os.path.join(dir_path, "..", "Demographics", "Outbreak_Villages")
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

# Parameter values
Taniuchi_village_files = search_files(
    r'/n/home04/weswong/multiscale_polio_evo/MultiscaleModeling/PopSim/Demographics/Matlab_Villages', '.pkl') #replace with location of files where the pregenerated Matlab village files are
idx = int(sys.argv[1])  # ranges from 1-20
strain_type = 'S2'

file = Taniuchi_village_files[idx]
print(Taniuchi_village_files)
print(file)
transmission_configs = {"fecal_oral_dose": 0,
                        "hh_transmission": 1,
                        "beta_hh": 0,
                        "bari_transmission": 1,
                        "beta_bari": 0,
                        "village_transmission": 1,
                        "beta_village": 0,
                        "age_assortivity": 0,
                        "beta_inter_village": 0,
                        "global_transmission": 0}

config.params.overwrite({"strain_type": strain_type})  # strain type of birthed individuals
config.params.overwrite(transmission_configs)

with open(file, 'rb') as fp:
    D = dill.load(fp)
write_dill(D, basename='Outbreak_Villages_0', output_dir=output_dir)

individuals = []
hh_ids = []
for V in D.villages:
    for bari in V.baris.values():
        for hhid in bari.households:
            hh_ids.append(hhid)
    individuals += V.return_individuals()
start_iid = np.max([i.id for i in individuals])
start_hid = np.max(hh_ids)
Individual.override_id_generator(start_iid + 10)
Household.override_hid_generator(start_hid + 10)

for year in range(1, 41):
    print(year)
    for V in D.villages:
        V.update(1)
    if year in [1, 2, 3, 4, 5, 10, 20, 40]:
        write_dill(D, basename='Outbreak_Villages_{idx}_{year}'.format(year=year, idx=str(idx)), output_dir=output_dir)
