import numpy as np
import itertools
import sys
import json
import shutil,os,glob
from collections import defaultdict
import dill
import Utility

collection = {} #temporary collection container
prevalence_dict = defaultdict(list)
r0 = defaultdict(lambda: defaultdict(list))
sim_cohorts = defaultdict(lambda: defaultdict(list))
n_transmission_levels= defaultdict(list)
for directory in glob.glob('*/'):
    print(directory)
    config = json.load(open(directory + 'config.json'))['config'][0]
    param_string = '{beta_hh}_{beta_bari}_{beta_village}_{beta_inter_village}_{dose}_{age_assortivity}_{recomb}_{genome_size}'.format(beta_hh = config['beta_hh'],
                                                                                                                beta_bari = config['beta_bari'],
                                                                                                                beta_village = config['beta_village'],
                                                                                                                beta_inter_village = config['beta_inter_village'],
                                                                                                                dose = config['fecal_oral_dose'],
                                                                                                                age_assortivity = config['age_assortivity'],
                                                                                                                recomb = config['recomb_rate'], genome_size = config['genome_size'])
    
    if param_string not in collection:
        collection[param_string] = defaultdict(list)
    
    file = directory + 'genetics.json'
    genetics_dict = json.load(open(file))
    

    for key in genetics_dict:
        collection[param_string][key] += genetics_dict[key]

for param_string in collection:
    destination_filename = 'genetics_{param_string}_prevalences.json'.format(param_string=param_string)
        
    with open(destination_filename, 'w') as fp:
        json.dump(collection[param_string], fp)
        
