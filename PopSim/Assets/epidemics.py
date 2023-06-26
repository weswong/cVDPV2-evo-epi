import os
import sys
from collections import defaultdict
from Demographics import Demographics, write_dill #this is necessary
import dill
import config
import json
import numpy as np

dir_path = os.path.dirname(os.path.realpath(__file__))
output_dir = os.path.join(dir_path, "..", "param_search") #sweep of fecal oral
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

def search_files(directory='.', extension='.pkl'):
    found_files  = []
    extension = extension.lower()
    for dirpath, dirnames, files in os.walk(directory):
        for name in files:
            if extension and name.lower().endswith(extension):
                found_files.append(os.path.join(dirpath, name))
       
    return found_files

idx = int(sys.argv[1]) - 1

#this is a set of demographics files that are generated prior to running. See
Matlab_village_files = search_files(r'/n/home04/weswong/multiscale_polio_evo/MultiscaleModeling/PopSim/Demographics/Matlab_Villages', '.pkl')

fecal_oral_doses = []
base_ranges = [_/10. for _ in range(10, 100, 5)]
dose_range = [1e-8, 1e-7]
for y in dose_range:
    for x in base_ranges:
        fecal_oral_doses.append(x*y)
fecal_oral_dose = fecal_oral_doses[idx]
transmission_configs = {"fecal_oral_dose" : fecal_oral_dose,
                        "hh_transmission" : 1,
                    "beta_hh" : 1,
                    "bari_transmission" : 1,
                    "beta_bari" : 15,
                    "village_transmission" : 1,
                    "beta_village" : 4,
                    "age_assortivity": 0,
                    "beta_inter_village":2,
                    "global_transmission" : 0}

config.params.overwrite(transmission_configs)
prevalences = defaultdict(list)
r0, transmissions = [], []
sim_cohort_individuals = defaultdict(list)
sim_cohort_incidences = defaultdict(list)
transmission_levels = defaultdict(list)
genetics = defaultdict(list)

for _ in range(50):#(30):
    random_idx = np.random.randint(0, len(Matlab_village_files) -1 )
    village_object_file = Matlab_village_files[random_idx]
    cohort_individuals = defaultdict(lambda: 0)
    
    
    with open(village_object_file, 'rb') as fp:
        D = dill.load(fp)

    D.Taniuchi_study(config.params)
    for key in D.sim_prevalences:
        prevalences[key].append(D.sim_prevalences[key])
        
    for key in D.cohort_incidences:
        sim_cohort_incidences[key].append(D.cohort_incidences[key])
    for gene in D.sim_genetics:
        genetics[gene].append(D.sim_genetics[gene])

    for V in D.villages:
        for key in V.infants:
            cohort_individuals['infants_' + key] += len([i.id for i in V.infants[key]])
        for key in V.hh_contacts:
            cohort_individuals['hh_' + key] += len([i.id for i in V.hh_contacts[key]])

    for key in D.n_events_per_level:
        transmission_levels[key].append(D.n_events_per_level[key])

    for key in cohort_individuals:
        sim_cohort_individuals[key].append(cohort_individuals[key])

fout = '{f}'.format(f = fecal_oral_dose)
with open(os.path.join(output_dir, '') + '{x}.json'.format(x=fout), 'w') as fp:
    json.dump([prevalences, sim_cohort_individuals, sim_cohort_incidences, transmission_levels], fp)
with open(os.path.join(output_dir, '') + '{x}_genetics.json'.format(x=fout), 'w') as fp:
    json.dump(genetics, fp)