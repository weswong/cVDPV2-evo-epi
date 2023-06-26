import numpy as np
import itertools
import sys
import json
import shutil,os,glob
from collections import defaultdict
import dill
import json
from collections import defaultdict


import Utility

results = defaultdict(lambda: defaultdict(list))
genetics_results = defaultdict(lambda: defaultdict(list))
incidence_results = defaultdict(lambda: defaultdict(list))
immunity_results = defaultdict(lambda: defaultdict(list))
for file in glob.glob('(')
        
    file = directory + 'epidemic.json'
    file_data = json.load(open(file, 'rb'))
    file_data = json.loads(file_data)
    outbreak_data = results[param_string]
    for key in file_data:
        if key != 'shed_durations':
            for repetition in file_data[key]:
                    outbreak_data[key].append([float(x) for x in repetition])
        else:
            outbreak_data[key] += file_data['shed_durations']
                
                
    genetics_file = directory + 'genetics.json'
    genetics_data = json.load(open(genetics_file))
    outbreak_genetics = genetics_results[param_string]
        
    for key in genetics_data:
        for repetition in genetics_data[key]:
            outbreak_genetics[key].append([x for x in repetition])
                
    immunity_file = directory + 'immunity.json'
    immunity_data = dill.load(open(immunity_file, 'rb'))
    outbreak_immunity = immunity_results[param_string]
        
    for key in immunity_data:
        for repetition in immunity_data[key]:
            outbreak_immunity[key].append([float(x) for x in repetition])
            outbreak_immunity[key].append([float(x) for x in repetition])

for param_string in results:            
    with open('epidemic_{param_string}.json'.format(param_string=param_string), 'w') as fp:
        json.dump(results[param_string], fp)
    with open('genetics_{param_string}.json'.format(param_string=param_string), 'w') as fp:
        json.dump(genetics_results[param_string], fp)
    #with open('incidences_{param_string}.json'.format(param_string=param_string), 'w') as fp:
    #    json.dump(incidence_results[param_string], fp)
    with open('immunity_{param_string}.json'.format(param_string=param_string), 'w') as fp:
        json.dump(immunity_results[param_string], fp)
        
