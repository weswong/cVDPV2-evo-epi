import json
import sys

dose_multiplier = float(sys.argv[1])

configs = []
years = [1,2,3,4,5]
for year in years:
    for p in [x/10. for x in range(9)]:
        for strain_type in ['S2']:#, 'C1', 'C2']:
            if p == 0:
                initial_conditions = 0
            else:
                initial_conditions = 1
            for iteration in range(100): #50
                c = {'year': year,
                    'p' : p,
                    'initial_conditions': initial_conditions,
                    'evo_strain_type':ls
                     strain_type,
                    "fecal_oral_dose" : 4e-7 * dose_multiplier,
                    'iteration': iteration}
                configs.append(c)
print(len(configs))
with open('outbreak_configs_{dose_multiplier}.json'.format(dose_multiplier = dose_multiplier), 'w') as f:
    json.dump(configs, f)
