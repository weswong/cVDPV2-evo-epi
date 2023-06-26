import numpy as np
from collections import defaultdict, Counter
from Village import Village
import dill
from Utility import modified_counts, draw_potential_infected
import config
import sys
from Infection import ImmunoInfection

class Demographics:
    
    def __init__(self):
        self.set_stats()
    
    def set_stats(self, params= config.params):
        self.years = params.sim_expansion_years + params.sim_contraction_years
        self.sim_expansion_years = params.sim_expansion_years
        self.sim_contraction_years = params.sim_contraction_years
        self.n_initial_villages = params.n_initial_villages
        #self.marriage_differentials = defaultdict(list)
        self.n_individuals = np.zeros((self.n_initial_villages, self.years))
        self.n_infants = np.zeros((self.n_initial_villages, self.years))
        self.n_baris = np.zeros((self.n_initial_villages, self.years))
        self.n_households = np.zeros((self.n_initial_villages, self.years))
        self.m_age_marriage = np.zeros((self.n_initial_villages, self.years))
        self.f_age_marriage = np.zeros((self.n_initial_villages, self.years))
        self.hh_sizes = defaultdict(lambda: np.zeros((self.n_initial_villages, 10)))
        self.age_pyramid = defaultdict(lambda: np.zeros((self.n_initial_villages, 20)))
        self.bari_hh = defaultdict(lambda: np.zeros((self.n_initial_villages, 10)))
        self.bari_n_individuals = defaultdict(lambda: np.zeros((self.n_initial_villages, 40)))
        self.villages = []
        self.p_no_infant_baris = defaultdict(list)
        self.infant_bari_count = defaultdict(lambda: np.zeros((self.n_initial_villages, 10)))
        self.community_bari_count = defaultdict(lambda: np.zeros((self.n_initial_villages, 10)))
        self.transmissions = defaultdict(list)
        self.sim_prevalences = defaultdict(list)
        self.historical_fertility_multiplier = params.historical_fertility_multiplier
     
    def fill_stats(self, V, iteration, year):
        self.p_no_infant_baris[iteration].append(V.noinfant_bari_prop())
        n_individuals, n_infants = V.n_individuals(infants = True)    
        self.n_individuals[iteration][year] = n_individuals
        self.n_infants[iteration][year] = n_infants
        self.n_households[iteration][year] = V.n_households()
        self.n_baris[iteration][year] = V.n_baris()
                
        individual_ages = []
        age_bins = range(5, 100, 5)
        for B in V.baris.values():
            for hh in B.households.values():
                individual_ages += [i.age for i in hh.individuals.values()]
        individual_age_counts = Counter(np.digitize(individual_ages, age_bins))
        modified_individual_age_counts = {}
        for _ in range(len(age_bins)):
            self.age_pyramid[year][iteration][_] = individual_age_counts.pop(_,0)
                    
                    
            hh_size_count = Counter(V.hh_sizes())
            self.hh_sizes[year][iteration] = modified_counts(hh_size_count)
                
            bari_n_hh_count = Counter(V.bari_n_households())
            self.bari_hh[year][iteration] = modified_counts(bari_n_hh_count)
                
            bari_nindividuals_count = Counter(V.bari_n_individuals())
            self.bari_n_individuals[year][iteration] = modified_counts(bari_nindividuals_count, 39)
                
                
            infant_bari_count = Counter(V.infant_bari_count())
            self.infant_bari_count[year][iteration] = modified_counts(infant_bari_count)    
    
    
    def initialize_default_villages(self, params = config.params): #initiate simulation such that all villages created have the same initial conditions
        params.overwrite({'historical_fertility_scale_fn': lambda x: params.historical_fertility_multiplier*(0.085*x + 0.18),
                              'm_death_age_sampler': params.historical_m_death_age_sampler,
                              'f_death_age_sampler': params.historical_f_death_age_sampler})
        print('Creating Villages')
        for iteration in range(self.n_initial_villages):
            V= Village.initialize_village(params.n_initial_baris, params.n_initial_hh_per_bari)
            self.villages.append(V)
        year =0
        
        print('Simulating Historical Setting')
        for _ in range(self.sim_expansion_years):
            print('Historical Period, Year ' + str(_))
            for village_idx,V in enumerate(self.villages):
                V.update(1)
                self.fill_stats(V, village_idx, year)
            year += 1

        print('Simulating demographic transition')
        for t in range(self.sim_contraction_years):
            print('Transition Year ' + str(t))
            scale = 1/self.sim_contraction_years*(t+1)
            params.overwrite({'historical_fertility_scale_fn' : lambda x: scale*(1-params.historical_fertility_multiplier*(0.085*x + 0.18)) + params.historical_fertility_multiplier*(0.085*x + 0.18),
                                  'f_death_age_sampler': lambda x: scale * (params.f_death_age_sampler_2015(x) - params.historical_f_death_age_sampler(x)) + params.historical_f_death_age_sampler(x),
                                  'm_death_age_sampler': lambda x: scale*(params.m_death_age_sampler_2015(x) - params.historical_m_death_age_sampler(x)) + params.historical_m_death_age_sampler(x)
                            })
            for village_idx, V in enumerate(self.villages):
                V.update(1)
                self.fill_stats(V,village_idx, year)
            year += 1
        for village_idx, V in enumerate(self.villages):
            print('Vaccinating Village {n}'.format(n = village_idx))
            V.vaccination_campaign(config.params.vaccination_campaign)
        
        
    def initialize_Taniuchi_villages(self, start_conditions_json='Taniuchi_village_start.json', params = config.params,
                                     demographic_decline = True, initialize_pre_existing_immunity= True): #initiate simulation to represent the villages used in the Taniuchi study
        params.set_start_conditions(start_conditions_json)
        self.set_stats(params)
        start_conditions = params.start_conditions #a dictionary where {target_village_size: [n_initial_baris, years, iterations]} assumes 1 household per bari
        params.overwrite({'historical_fertility_scale_fn': lambda x: params.historical_fertility_multiplier*(0.085*x + 0.18),
                              'm_death_age_sampler': params.historical_m_death_age_sampler,
                              'f_death_age_sampler': params.historical_f_death_age_sampler})
        print('Initializing Villages')
        for target_village_size in start_conditions:
            n_initial_baris, repetitions = int(start_conditions[target_village_size][0]),int(start_conditions[target_village_size][1])
            for rep in range(repetitions):
                V = Village.initialize_village(n_initial_baris, 1)
                self.villages.append(V)
            
        print('Simulating Historical Setting')
        year =0 
        for _ in range(self.sim_expansion_years):
            for village_idx,V in enumerate(self.villages):
                V.update(1)
                self.fill_stats(V, village_idx, year)
            year += 1
        
        if demographic_decline == True:
            print('Simulating demographic transition')
            for t in range(self.sim_contraction_years):
                scale = 1/self.sim_contraction_years*(t+1)
                params.overwrite({'historical_fertility_scale_fn' : lambda x: scale*(1-params.historical_fertility_multiplier*(0.085*x + 0.18)) + params.historical_fertility_multiplier*(0.085*x + 0.18),
                                      'f_death_age_sampler': lambda x: scale * (params.f_death_age_sampler_2015(x) - params.historical_f_death_age_sampler(x)) + params.historical_f_death_age_sampler(x),
                                      'm_death_age_sampler': lambda x: scale*(params.m_death_age_sampler_2015(x) - params.historical_m_death_age_sampler(x)) + params.historical_m_death_age_sampler(x)
                                })
                for village_idx, V in enumerate(self.villages):
                    V.update(1)
                    self.fill_stats(V,village_idx, year)
                year += 1
        else:
            for t in range(self.sim_contraction_years):
                for village_idx, V in enumerate(self.villages):
                    V.update(1)
                    self.fill_stats(V,village_idx, year)
                year += 1
        
        print('Vaccinating Villages')
        for V in self.villages:
            V.vaccination_campaign(params.vaccination_campaign, initialize_pre_existing_immunity = initialize_pre_existing_immunity)
        
            
            
        params.overwrite({'historical_fertility_scale_fn': lambda x: 1,
                              'm_death_age_sampler': params.m_death_age_sampler_2015,
                              'f_death_age_sampler': params.f_death_age_sampler_2015})

    
        
    
    def inter_village_transmission(self, transmissions, params = config.params):
        transmissions_dict = defaultdict(list)
        for Tx in transmissions:
            transmissions_dict[Tx.village_id].append(Tx)
        if params.beta_inter_village > 0:
            for village_id in transmissions_dict:
                #print(village_id, len(transmissions_dict[village_id]))
                source_village = village_id
                mask = self.village_tracker != source_village
                infectees = self.individuals[mask]
                infectee_ages = self.individual_ages[mask]

                ego_ages, viral_doses, tx_tracker, genome_tracker = [], [], [],[]
                for idx, Tx in enumerate(transmissions_dict[village_id]):      
                    ego_ages.append(Tx.source_age)
                    viral_doses += [Tx.dose * params.fecal_oral_dose] * params.beta_inter_village #python list multiplication does not multiply each element but replicates it
                    genome_tracker += [Tx.genome] * params.beta_inter_village
                    tx_tracker += [(village_id, idx)] * params.beta_inter_village
                ego_age_counts = Counter(ego_ages)

                if len(infectees) > 0:
                    infectee_idxes = []

                    for ego_age in ego_age_counts: #done this way to prevent mass look up of the age assortivity matrix....
                        number_of_samples = ego_age_counts[ego_age] * params.beta_inter_village
                        infectee_idxes += list(draw_potential_infected(infectees, params.age_assortivity_non_uniform,
                                    number_of_samples, params.village_age_assortivity_matrix[ego_age][infectee_ages]))

                    #print('idx viral dose counts', len(infectee_idxes), len(viral_doses)) 
                    #print('idx viral dose counts', len(infectee_idxes), len(viral_doses))

                    for infectee_idx, dose, tx_index, genome in zip(infectee_idxes, viral_doses, tx_tracker, genome_tracker):
                        infected = infectees[infectee_idx]
                        infected_status, prechallenge_immunity = infected.receive_transmission_dose(dose, genome, params = params)
                        if infected_status:
                            transmissions_dict[tx_index[0]][tx_index[1]].secondary_cases['total'] += 1
                            if infected.age <= 5:
                                transmissions_dict[tx_index[0]][tx_index[1]].secondary_cases['0-5'] += 1
                            elif infected.age <= 15:
                                transmissions_dict[tx_index[0]][tx_index[1]].secondary_cases['5-15'] += 1
                            else:
                                transmissions_dict[tx_index[0]][tx_index[1]].secondary_cases['>15'] += 1
                            transmissions_dict[tx_index[0]][tx_index[1]].tx_history['inter_village'].append((infected.id, infected.age, prechallenge_immunity, infected.times_infected))



    def global_transmission(self, transmissions, params = config.params):
        viral_doses = []
        potential_infectee_indexes = range(len(self.individuals))
        ego_ages, viral_doses, tx_tracker, genome_tracker, meta_tracker = [], [], [],[], []
        for idx, Tx in enumerate(transmissions):
            ego_ages.append(Tx.source_age)
            viral_doses += [Tx.dose * params.fecal_oral_dose] * params.beta_global
            genome_tracker += [Tx.genome] * params.beta_global
            tx_tracker += [idx] * params.beta_global
            meta_tracker += [(Tx.hhid, Tx.bari_id, Tx.village_id)]* params.beta_global
        ego_age_counts = Counter(ego_ages)
        infectee_idxes = []
        
        for ego_age in ego_age_counts:
            number_of_samples = ego_age_counts[ego_age] * params.beta_global
            infectee_idxes+= list(draw_potential_infected(potential_infectee_indexes, params.age_assortivity_non_uniform,
                                                           number_of_samples, params.village_age_assortivity_matrix[ego_age][self.individual_ages]))
        
        for infectee_idx, dose, tx_idx, genome in zip(infectee_idxes, viral_doses, tx_tracker, genome_tracker):
            infected_status, prechallenge_immunity = self.individuals[infectee_idx].receive_transmission_dose(dose, genome, params = params)
            #print(self.individuals[infectee_idx].age, self.individuals[infectee_idx].infection.current_immunity, self.individuals[infectee_idx].infection.p_infection(dose))
            
            if infected_status:
                transmissions[tx_idx].secondary_cases['total'] += 1
                if self.individuals[infectee_idx].age <= 5:
                    transmissions[tx_idx].secondary_cases['0-5'] += 1
                elif self.individuals[infectee_idx].age  <= 15:
                    transmissions[tx_idx].secondary_cases['5-15'] += 1
                else:
                    transmissions[tx_idx].secondary_cases['>15'] += 1
                #print(self.individuals[infectee_idx].id, self.village_tracker[infectee_idx], self.individuals[infectee_idx].isalive, self.individuals[infectee_idx].age)
                infectee_hh, infectee_bari, infectee_v = self.track_individuals(self.individuals[infectee_idx].id)
                source_hh, source_bari, source_v = meta_tracker[tx_idx]
                #print(source_hh, source_bari, source_v)
                #print(infectee_hh, infectee_bari, infectee_v)
                
                if int(infectee_v) == int(source_v):
                    if int(infectee_bari) == int(source_bari):
                        if int(infectee_hh) == int(source_hh):
                            level = 'hh'
                        else:
                            level = 'bari'
                    else:
                        level = 'village'
                else:
                    level = 'inter_village'
                #print(level)
                    
                transmissions[tx_idx].tx_history[level].append((self.individuals[infectee_idx].id, self.individuals[infectee_idx].age,prechallenge_immunity, self.individuals[infectee_idx].times_infected))
        
    def track_individuals(self,iid):
        for V in self.villages:
            if iid in [i.id for i in V.return_individuals()]:
                for bari in V.baris.values():
                    if iid in [i.id for i in bari.return_individuals()]:
                        for hh in bari.households.values():
                            if iid in [i.id for i in hh.individuals.values()]:
                                return [hh.id, bari.id, V.id]
            
    def find_individuals(self):
        self.village_tracker = []
        all_bari_individuals = []
        for V in self.villages:
            for bari in V.baris.values():
                bari_individuals = bari.return_individuals()
                all_bari_individuals += bari_individuals
                self.village_tracker += len(bari_individuals) * [V.id]

        self.individual_ages = np.asarray([int(round(I.age, 0)) if I.age < 80 else 80 for I in all_bari_individuals])
        self.individuals = np.asarray(all_bari_individuals)
        self.village_tracker = np.asarray(self.village_tracker)
    
    def set_outbreak_conditions(self, strain_type = 'S2', params = config.params):
        individuals = []
        for V in self.villages:
            individuals += V.return_individuals()
        
        for i in individuals:
            if i.age < 1:
                i.infection = ImmunoInfection.initialize_infant_immunity('bOPV', params = params)
            else:
                i.infection = ImmunoInfection.initialize_noninfant_immunity(i.age, 'bOPV', params = params) 
            i.infection.strain_type = strain_type
                    
    
    
    def outbreak(self, p = 0.0, children = False, single_infant = False, params = config.params, redo_vaccination = False, strain_type = 'S2'):
        self.find_individuals()    #find all the individuals in the simulation
        self.village_counts = defaultdict(list)
        self.incidence = defaultdict(list)
        self.cohort_incidence = defaultdict(list)
        r0_tracker = defaultdict(list)
        self.transmission_tracker = {}
        self.infection_ids = defaultdict(list)
        self.cohort_incidences = defaultdict(lambda: np.zeros(161))
        self.cumulative_cohort_incidences = defaultdict(lambda: np.zeros(161))
        self.cohort_individuals = defaultdict(list)
        self.sim_prevalences = defaultdict(list)
        self.sim_genetics = defaultdict(list)
        self.sim_immunity = defaultdict(list)
        
        sim_infants = []
        prevalence, transmissions, demographics, genetics = defaultdict(list),[],[], defaultdict(list)
        immunity = defaultdict(list)
        
        individuals = []
        for V in self.villages:
            individuals += V.return_individuals()
        
        if single_infant == True:
            infants = []
            for i in individuals:
                if i.age <=1:
                    infants.append(i)
                i.times_infected = 0
                
            n_targets = len(infants)
            random_idx = np.random.choice([x for x in range(n_targets)], 1)[0]
            random_index_case = infants[random_idx]
            random_index_case.infection = ImmunoInfection.initialize_null_immunity(strain_type)
            random_index_case.infection.attempt_infection(params=params)#genome= {'A481G': 1, 'U2909C': 1, 'U398C': 1, 'nonsyn': 9, 'nonsyn_neutral': 0, 'syn' : 0}, params = params)
            idx_shed_duration = random_index_case.infection.shed_duration
            
        elif children == True:
            children = []
            for i in individuals:
                if i.age <=5:
                    children.append(i)
                i.times_infected = 0
                
            n_children = len(children)
            n_indexes = np.random.choice([x for x in range(n_children)], int(round(p * n_children, 0))) 
            for idx in n_indexes:
                random_index_case = children[idx]
                random_index_case.infection = ImmunoInfection.initialize_null_immunity(strain_type)
                random_index_case.infection.attempt_infection(params = params)
        
        for i in individuals:
            prevalence['total'].append(i.is_transmitting)
            immunity['low'].append(i.infection.current_immunity <= 8)
            immunity['mid'].append((i.infection.current_immunity > 8) and (i.infection.current_immunity < 256))
            immunity['high'].append(i.infection.current_immunity >= 256)
            if i.age <= 5:
                prevalence['0-5'].append(i.is_transmitting)
                immunity['low_0-5'].append(i.infection.current_immunity <= 8)
                immunity['mid_0-5'].append((i.infection.current_immunity > 8) and (i.infection.current_immunity < 256))
                immunity['high_0-5'].append(i.infection.current_immunity >= 256)
            elif i.age <= 15:
                prevalence['5-15'].append(i.is_transmitting)
                immunity['low_5-15'].append(i.infection.current_immunity <= 8)
                immunity['mid_5-15'].append((i.infection.current_immunity > 8) and (i.infection.current_immunity < 256))
                immunity['high_5-15'].append(i.infection.current_immunity >= 256)
            else:
                prevalence['>15'].append(i.is_transmitting)
                immunity['low_15+'].append(i.infection.current_immunity <= 8)
                immunity['mid_15+'].append((i.infection.current_immunity > 8) and (i.infection.current_immunity < 256))
                immunity['high_15+'].append(i.infection.current_immunity >= 256)
            
                    
            if i.is_transmitting:
                for mutation in params.loci:
                    genetics[mutation].append(i.infection.genome[mutation])
                genetics['transmission_fitness'].append(i.infection.calculate_transmission_fitness_self(params))
                genetics['duration_fitness'].append(i.infection.naiive_shed_duration_reporter)
                genetics['vp1'].append(i.infection.calculate_vp1_mutations())

                        
                        
        secondary_cases=defaultdict(lambda: 0)
        for tx in transmissions:
            for key in tx.secondary_cases:
                secondary_cases[key] += tx.secondary_cases[key]
        for key in ['total', '0-5', '5-15','>15']:
            self.incidence[key].append(secondary_cases[key])
        
        prevalence['n_total'].append(len(prevalence['total']))
        prevalence['n_0-5'].append(len(prevalence['0-5']))
        prevalence['n_5-15'].append(len(prevalence['5-15']))
        prevalence['n_>15'].append(len(prevalence['>15']))
        for key in prevalence:
            self.sim_prevalences[key].append(np.sum(prevalence[key]))
            print(key, np.sum(prevalence[key]))
        for gene in genetics:
            self.sim_genetics[gene].append([np.mean(genetics[gene]), np.percentile(genetics[gene], 5), np.percentile(genetics[gene], 95)])
        for key in immunity:
            self.sim_immunity[key].append(np.mean(immunity[key]))
        
            
        for day in range(365 * 3):
            prevalence, transmissions, demographics, genetics = defaultdict(list),[],[], defaultdict(list)
            immunity = defaultdict(list)
            print('Day {day}: Updating Villages'.format(day = day))

            # !!!!!!!!! DEBUG !!!!!!!!!!!!!!!!
            # if day == 10:
            #     sys.exit(1)
            # !!!!!!!!!

            for V in self.villages:
                V.update(1/365)
                transmissions += V.transmissions
                demographics += V.demographics_tracker
            print('Index case is still shedding: ' + str(random_index_case.infection.is_shed))
               
            print('Day {day}: Updating Matlab individual stats'.format(day = day))
            self.find_individuals()
                
            if params.global_transmission == 0:
                self.inter_village_transmission(transmissions, params)
                print('Day {day}: Initiating Inter-Village Transmission'.format(day = day) )
            elif params.global_transmission == 1:
                print('global')
                self.global_transmission(transmissions, params)
                
            else:
                continue #if global_transmission == 2, it is handled at the village level for "island model"
                
            print('Day {n}'.format(n = day))
            for i in self.individuals:
                prevalence['total'].append(i.is_transmitting)
                immunity['low'].append(i.infection.current_immunity <= 8)
                immunity['mid'].append((i.infection.current_immunity > 8) and (i.infection.current_immunity < 256))
                immunity['high'].append(i.infection.current_immunity >= 256)
                if i.age <= 5:
                    prevalence['0-5'].append(i.is_transmitting)
                    immunity['low_0-5'].append(i.infection.current_immunity <= 8)
                    immunity['mid_0-5'].append((i.infection.current_immunity > 8) and (i.infection.current_immunity < 256))
                    immunity['high_0-5'].append(i.infection.current_immunity >= 256)
                elif i.age <= 15:
                    prevalence['5-15'].append(i.is_transmitting)
                    immunity['low_5-15'].append(i.infection.current_immunity <= 8)
                    immunity['mid_5-15'].append((i.infection.current_immunity > 8) and (i.infection.current_immunity < 256))
                    immunity['high_5-15'].append(i.infection.current_immunity >= 256)
                else:
                    prevalence['>15'].append(i.is_transmitting)
                    immunity['low_15+'].append(i.infection.current_immunity <= 8)
                    immunity['mid_15+'].append((i.infection.current_immunity > 8) and (i.infection.current_immunity < 256))
                    immunity['high_15+'].append(i.infection.current_immunity >= 256)
                    
                immunity['low'].append(i.infection.current_immunity <= 8)
                immunity['mid'].append((i.infection.current_immunity > 8) and (i.infection.current_immunity < 256))
                immunity['high'].append(i.infection.current_immunity >= 256)
                    
                
                if i.is_transmitting:
                    for mutation in params.loci:
                        genetics[mutation].append(i.infection.genome[mutation])
                    genetics['transmission_fitness'].append(i.infection.calculate_transmission_fitness_self(params))
                    genetics['duration_fitness'].append(i.infection.calculate_naiive_shed_duration(params))
                    genetics['vp1'].append(i.infection.calculate_vp1_mutations())
                        
                   
            
            for gene in genetics:
                self.sim_genetics[gene].append([np.mean(genetics[gene]), np.percentile(genetics[gene], 5), np.percentile(genetics[gene], 95)])
            
            secondary_cases=defaultdict(lambda: 0)
            for tx in transmissions:
                for key in tx.secondary_cases:
                    secondary_cases[key] += tx.secondary_cases[key]
            for key in ['total', '0-5', '5-15','>15']:
                self.incidence[key].append(secondary_cases[key])
                
            for key in immunity:
                self.sim_immunity[key].append(np.mean(immunity[key]))
                    
            
            self.transmission_tracker[day] = transmissions
            prevalence['n_total'].append(len(prevalence['total']))
            prevalence['n_0-5'].append(len(prevalence['0-5']))
            prevalence['n_5-15'].append(len(prevalence['5-15']))
            prevalence['n_>15'].append(len(prevalence['>15']))
            for key in prevalence:
                self.sim_prevalences[key].append(np.sum(prevalence[key]))
                print(key, np.sum(prevalence[key]))
            
            
            if np.sum(prevalence['total']) ==0:
                break
            
        n_events_per_level, tx_age_matrix   = extract_transmissions(self.transmission_tracker, final_day = day)
        self.n_events_per_level = n_events_per_level
        
        if single_infant:
            return idx_shed_duration
            
    def Taniuchi_study(self, params = config.params, redo_vaccination = False):
        if redo_vaccination:
            for V in self.villages:
                V.vaccination_campaign(params.vaccination_campaign, initialize_pre_existing_immunity= True)
                print(V.infants)
                print(V.hh_contacts)
                print('Vaccination of Village {n} Complete'.format(n= V.id))
        self.find_individuals()    #find all the individuals in the simulation
        self.village_counts = defaultdict(list)
        self.incidence = defaultdict(list)
        self.cohort_incidence = defaultdict(list)
        r0_tracker = defaultdict(list)
        self.transmission_tracker = {}
        self.infection_ids = defaultdict(list)
        self.cohort_incidences = defaultdict(lambda: np.zeros(161))
        self.cumulative_cohort_incidences = defaultdict(lambda: np.zeros(161))
        self.cohort_individuals = defaultdict(list)
        self.sim_prevalences = defaultdict(list)
        self.sim_genetics = defaultdict(list) # for A481G, U2909C, U398C, nonsyn
        
        for V in self.villages:
            for individual in V.return_individuals():
                if individual.infection.is_shed:
                    individual.n_infection = 1
                    print(individual.id, individual.infection.genome)
                else:
                    individual.n_infection = 0
        
        for day in range(23*7):
            prevalence, transmissions, demographics, genetics = defaultdict(list),[],[], defaultdict(list)
            print('Day {day}: Updating Villages'.format(day = day))
            for V in self.villages:
                V.update(1/365)
                transmissions += V.transmissions
                demographics += V.demographics_tracker
               
            print('Day {day}: Updating Matlab individual stats'.format(day = day))
            self.find_individuals()
                
            print('Day {day}: Initiating Inter-Village Transmission'.format(day = day) )
            if params.global_transmission == 0:
                self.inter_village_transmission(transmissions, params)
            elif params.global_transmission == 1:
                self.global_transmission(transmissions, params)
            else:
                continue #if global_transmission == 2, it is handled at the village level for "island model"

                
            print('Day {day}: Recording Stats'.format(day = day))  
            if day in 7*np.asarray([0,1,2,3,4,5,6,7,8,9,10, 14, 18, 22]):
                print('Day {n}'.format(n = day))
                for V in self.villages:
                    surveillance_cohort = []
                    for key in V.infants:
                        prevalence['infants_' + key] += [i.is_transmitting for i in V.infants[key]]
                        self.infection_ids['infants_' + key].append([i.id for i in V.infants[key]])
                        
                        surveillance_cohort += prevalence['infants_' + key]
                    for key in V.hh_contacts:
                        prevalence['hh_'+str(key)] += [i.is_transmitting for i in V.hh_contacts[key]]
                        surveillance_cohort += prevalence['hh_'+str(key)]
                        self.infection_ids['hh_' + key].append([i.id for i in V.hh_contacts[key]])
                    
                    village_prevalence_array= []
                    for i in V.return_individuals():
                        prevalence['total'].append(i.is_transmitting)
                        if i.is_transmitting:
                            for mutation in params.loci:
                                genetics[mutation].append(i.infection.genome[mutation])
                            genetics['transmission_fitness'].append(i.infection.calculate_transmission_fitness_self(params))
                            genetics['duration_fitness'].append(i.infection.calculate_naiive_shed_duration(params))
                            genetics['vp1'].append(i.infection.calculate_vp1_mutations())

                    self.village_counts[V.id].append([int(np.sum(surveillance_cohort)), len(surveillance_cohort), V.n_individuals()])
                
                
                
            for key in prevalence:
                self.sim_prevalences[key].append(np.mean(prevalence[key]))
            for gene in genetics:
                self.sim_genetics[gene].append([np.mean(genetics[gene]), np.percentile(genetics[gene], 5), np.percentile(genetics[gene], 95)])
                
            for tx in transmissions:
                self.incidence[day].append(tx.secondary_cases['total'])
                r0_tracker[tx.source_iid].append((day, tx.secondary_cases['total']))
            self.transmission_tracker[day] = transmissions
        

        for V in self.villages:
            for key in V.infants:
                self.cohort_individuals['infants_' + key] += [i.id for i in V.infants[key]]
            for key in V.hh_contacts:
                self.cohort_individuals['hh_' + key] += [i.id for i in V.hh_contacts[key]]

        for day in range(23*7):
            for Tx in [Tx for Tx in self.transmission_tracker[day] if Tx.secondary_cases != 0]:
                for cohort in self.cohort_individuals:
                    for level in Tx.tx_history:
                        for event in Tx.tx_history[level]:
                            infected_id = event[0]
                            if infected_id in self.cohort_individuals[cohort]:
                                self.cohort_incidences[cohort][day] += 1
        
        for cohort in self.cohort_individuals:
            self.cohort_incidences[cohort] = list(self.cohort_incidences[cohort])
        
        n_events_per_level, tx_age_matrix   = extract_transmissions(self.transmission_tracker)
        self.n_events_per_level = n_events_per_level
        
        self.r0 = []
        for infection_id in r0_tracker:
            t_0, incidence_count = r0_tracker[infection_id][0]
            for tracker in r0_tracker[infection_id][1:]:
                t_1, incidence = tracker
                if t_0 + 1 == t_1:
                    incidence_count += incidence
                else:
                    self.r0.append(incidence_count)
                    incidence_count = incidence
                t_0 = t_1
            self.r0.append(incidence_count)
        self.r0 = Counter(self.r0)
        
            
            
def extract_transmissions(transmissions, final_day = 161):
    age_bins = range(0,85,5)
    n_events_per_level = {}
    tx_age_matrix = defaultdict(lambda: defaultdict(lambda: np.zeros(18)))
    n_events = 0
    
    for day in range(final_day):
        for Tx in transmissions[day]:
            if Tx.secondary_cases != 0:
                digitized_source_age = int(np.digitize(Tx.source_age, age_bins))
                for level in Tx.tx_history:
                    for event in Tx.tx_history[level]:
                        n_events += 1
                        digitized_ego_age = int(np.digitize(event[1], age_bins))
                        tx_age_matrix[level][digitized_source_age][digitized_ego_age] += 1
    
    for level in tx_age_matrix:
        level_sum = np.sum(np.asarray(list(tx_age_matrix[level].values()))) #events per level
        n_events_per_level[level] = level_sum
        for digitized_ego_age in tx_age_matrix[level]:
            tx_age_matrix[level][digitized_ego_age] = tx_age_matrix[level][digitized_ego_age] / level_sum
    
    return n_events_per_level, tx_age_matrix        
        
def write_dill(object,  basename ='village_simulations', output_dir='../sims'):
        with open(output_dir + '/{basename}.pkl'.format(basename = basename), 'wb') as f:
            dill.dump(object, f)
        f.close()

            
if __name__=='__main__':
    import os
    import sys
    import json
    dir_path = os.path.dirname(os.path.realpath(__file__))
    D  = Demographics()
    #D.initialize_default_villages()
    #D.Taniuchi_study()
    #start_time = time.time()
    #D.initialize_Taniuchi_villages()
    D.initialize_Taniuchi_villages(dir_path + '/Taniuchi_village_start.json', initialize_pre_existing_immunity=True,
                                   demographic_decline=True)
    #end_time = time.time()
    #print("Elapsed time was %g seconds" % (end_time - start_time), flush = True)
    if len(sys.argv) > 1:
        file_index = str(sys.argv[1])
    else:
        file_index = str(0)
    write_dill(D, basename = 'Matlab_Village_{i}'.format(i=file_index), output_dir = '/n/home04/weswong/multiscale_polio_evo/MultiscaleModeling/PopSim/Demographics')