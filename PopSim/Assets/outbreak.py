import os
import sys
from collections import defaultdict
import dill
import config
import json
import numpy as np
from Individual import Individual
from Household import Household


    
def search_files(directory='.', extension='.pkl'):
    found_files  = []
    extension = extension.lower()
    for dirpath, dirnames, files in os.walk(directory):
        for name in files:
            if extension and name.lower().endswith(extension):
                found_files.append(os.path.join(dirpath, name))
    return found_files

param_idx = int(sys.argv[1])
iter2 = sys.argv[2]
dose_multiplier = float(sys.argv[3])

param_file = '/n/home04/weswong/multiscale_polio_evo/MultiscaleModeling/PopSim/outbreak_configs_v2_{dose_multiplier}.json'.format(dose_multiplier = dose_multiplier)

dir_path = os.path.dirname(os.path.realpath(__file__))
output_dir = os.path.join(dir_path, "..",  "output", 'dose_{d}'.format(d = dose_multiplier))
if not os.path.exists(output_dir):
    os.mkdir(output_dir)


with open(param_file) as fin:
    parameters = json.load(fin)

p = parameters[param_idx]['p']
strain_type = parameters[param_idx]['evo_strain_type']
iteration = parameters[param_idx]['iteration']
year = parameters[param_idx]['year']
#p, strain_type, year, iteration = parameters[year][param_idx]

fout = '{strain_type}_{p}_{year}_{iter}.{iter2}_{dose_multiplier}.'.format(strain_type = strain_type, p = str(p), year = str(year),
                                                                           iter = str(iteration), iter2 = str(iter2), dose_multiplier = str(dose_multiplier))
print(fout)
if float(p) == 0.0: #if p is zero, assume it is a point importation
    initial_conditions = 0
else:
    initial_conditions = 1 #otherwise initiate with p infections

Matlab_village_files = search_files(r'/n/home04/weswong/multiscale_polio_evo/MultiscaleModeling/PopSim/Demographics/Outbreak_Villages', '_{year}.pkl'.format(year = year))
print(Matlab_village_files)
transmission_configs = {"fecal_oral_dose" : 4e-7 * dose_multiplier,
                        "hh_transmission" : 1,
                    "beta_hh" : 1,
                    "bari_transmission" : 1,
                    "beta_bari" : 15,
                    "village_transmission" : 1,
                    "beta_village" : 4,
                    "age_assortivity": 0,
                    "beta_inter_village":2,
                    "global_transmission" : 0}


transmission_configs = transmission_configs
config.params.overwrite(transmission_configs)   

                              
prevalences = defaultdict(list)
genetics = defaultdict(list)
r0, transmissions = [], []
sim_cohort_individuals = defaultdict(list)
sim_cohort_incidences = defaultdict(list)
transmission_levels = defaultdict(list)
immunity = defaultdict(list)

for _ in range(1):
    print('Starting iteration: {i}'.format(i=_))
    random_idx = np.random.randint(0, len(Matlab_village_files) -1 )
    village_object_file = Matlab_village_files[random_idx] #Matlab_village_files[random_idx]
    cohort_individuals = defaultdict(list)
    with open(village_object_file, 'rb') as fp:
        D = dill.load(fp)
    individuals = []
    n_hh = 0
    hh_ids = []
    for V in D.villages:
        for bari in V.baris.values():
            for hhid in bari.households:
                hh_ids.append(hhid)
        individuals += V.return_individuals()
    start_iid = np.max([i.id for i in individuals])
    start_hid = np.max(hh_ids)
    Individual.override_id_generator(start_iid+10)
    Household.override_hid_generator(start_hid+10)

            


    if str(initial_conditions) == str(0):
        shed_duration = D.outbreak(p=0.0, children=False, single_infant=True, params=config.params, strain_type='S2')
    else:
        D.outbreak(p=p, children=True, single_infant=False, params=config.params, strain_type='S2')
        shed_duration = None

    for key in D.sim_prevalences:
        prevalences[key].append(D.sim_prevalences[key])
    for gene in D.sim_genetics:
        genetics[gene].append(D.sim_genetics[gene])
    for key in D.n_events_per_level:
        transmission_levels[key].append(D.n_events_per_level[key])
    for key in D.incidence:
        sim_cohort_incidences[key].append(D.incidence[key])
    for key in D.sim_immunity:
        immunity[key].append(D.sim_immunity[key])

class NumpyEncoder(json.JSONEncoder):
    """ Special json encoder for numpy types """

    def default(self, obj):
        if isinstance(obj, (np.int_, np.intc, np.intp, np.int8,
                            np.int16, np.int32, np.int64, np.uint8,
                            np.uint16, np.uint32, np.uint64)):
            return int(obj)
        elif isinstance(obj, (np.float_, np.float16, np.float32,
                              np.float64)):
            return float(obj)
        elif isinstance(obj, (np.ndarray,)):  #### This is the fix
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


with open(os.path.join(dir_path, "..", "output", 'dose_{d}_genetics'.format(d = dose_multiplier), '') + fout + '_genetics.json', 'w') as fp:
    json.dump(genetics, fp)

with open(os.path.join(dir_path, "..", "output", 'dose_{d}_genetics'.format(d = dose_multiplier), '') + fout +  '_epidemic.json', 'w') as fp:
    dumped = json.dumps(prevalences, cls=NumpyEncoder)
    json.dump(dumped, fp)

with open(os.path.join(dir_path, "..", "output", 'dose_{d}_genetics'.format(d = dose_multiplier), '') + fout +  '_immunity.json', 'w') as fp:
    dumped = json.dumps(immunity, cls=NumpyEncoder)
    json.dump(dumped, fp)
