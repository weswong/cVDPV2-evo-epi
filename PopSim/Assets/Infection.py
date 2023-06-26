import numpy as np
import math
import config
import json

class ImmunoInfection:
    '''typical infection/immune dynamics for poliovirus, this is also where pathogen-side genetics will be put (either that or split it into its own class)
    This is keyed towards polio -- currently completely avirulent!
    dt is measured in days!'''
    
    def __init__(self, strain_type):
        self.naiive = True
        self.IPV_naiive = False #received an IPV dose, but does not shed
        self.is_shed = False
        self.strain_type = strain_type
        self.t_infection = 0 #days, time since last infection
        self.current_immunity = 1 #Naiive
        self.prechallenge_immunity = 1
        self.postchallenge_peak_immunity = None
        self.shed_duration = 0
        self.genome = None
        self.optimum = False



    def initialize_evolution(self, genome = None, params = config.params): #genome = A481G, U2909C, U398C, NONSYN
        '''Use only for S2'''
        if not genome:
            self.genome = params.return_genome_template()
        else:
            self.genome = genome
        self.optimum = False
        #setting mutation wait times -- these mutations arise primarily through mutation pressure and are dfeined as within host fixations per bp per day
        self.mutation_wait_time = {}
        
        for mutation in self.genome:
            if mutation not in ['nonsyn', 'syn', 'cpg', 'nonsyn_neutral']:
                if self.genome[mutation] == 0:
                    self.mutation_wait_time[mutation] = self.mutate(params.fixation_rates[mutation])
                else:
                    self.mutation_wait_time[mutation] = math.inf
            else:
                self.mutation_wait_time[mutation] = self.mutate(params.fixation_rates[mutation])
        self.mutation_wait_time['recomb'] = self.mutate(params.recomb_rate)
        
        
        #self.check_epistasis()
        
        
    #def check_epistasis(self, params = config.params):
        #adaptive_mutations = [mutation for mutation in self.genome if mutation not in ['nonsyn', 'syn', 'cpg', 'nonsyn_neutral']]
        #if np.sum([self.genome[mutation] for mutation in adaptive_mutations]) == len(adaptive_mutations): #adaptive peak has been reached
        #    self.optimum = True
        #    self.mutation_wait_time['nonsyn_neutral'] = self.t_infection + self.mutate(params.fixation_rates['nonsyn_neutral'])
        
    def mutate(self, rate, dt = 1):
        if rate == 0:
            return 1e8
        else:
            return np.random.exponential(1/(dt * rate))
    
    def return_genome(self):
        return json.loads(json.dumps(self.genome))
        
        
    @classmethod
    def calculate_transmission_fitness(cls, beta, alleles, selection_coefficients):
        denominator = 1
        for a,s in zip(alleles, selection_coefficients):
            denominator *= (1 + s)**a
        w = beta / denominator        
        return w
    
    def calculate_transmission_fitness_self(self, params):
        '''to calculate the transmission fitness of the self infection'''
        beta = params.infectiousness_params['beta']
        selection_coefficients = [params.infectiousness_s[locus] for locus in params.loci]
        alleles = [self.genome[locus] for locus in params.loci]
        w = ImmunoInfection.calculate_transmission_fitness(beta, alleles, selection_coefficients)
        return w        
        
    def attempt_infection(self, dose=10**6, modifier = 1, genome = None, params = config.params): #genome = A481G, U2909C, U398C, NONSYN
        if not genome:
            genome = params.return_genome_template()
        loci = params.loci
        selection_coefficients = [params.infectiousness_s[locus] for locus in loci]
        alleles = [genome[locus] for locus in loci] 
        beta = params.infectiousness_params['beta']
        w = ImmunoInfection.calculate_transmission_fitness(beta, alleles, selection_coefficients)
        p_infection = (1 - (1 + dose/w)**(-params.infectiousness_params['alpha'] * (self.current_immunity)**-params.infectiousness_params['gamma'])) #* 0.554
    
        if (np.random.random() <  modifier * p_infection): #this actually causes the probability of infection to go beyond 1.0
            #however, np.random.random() only draws from 0-1 so it is the same as maxing it out
            self.is_shed=True
            self.naiive, self.IPV_naiive = False, False
            self.prechallenge_immunity = self.current_immunity
            self.t_infection = 0 #infection is reset
            self.calculate_peak_immunity()
            self.calculate_current_immunity()
            self.initialize_evolution(genome, params = params)
            self.calculate_shed_duration(params = params)

            return 1
        else:
            return 0
        
    def calculate_vp1_mutations(self):
        return 903/7440 * (self.genome['syn'] + self.genome['nonsyn'] + self.genome['nonsyn_neutral']) + self.genome['U2909C']
    
    
    def calculate_peak_immunity(self):
        '''immunity immediately post infection'''
        self.postchallenge_peak_immunity = self.prechallenge_immunity * max(1, self.theta_Nab()) #prevent immunity from decreasing due to challenge
    
    def calculate_current_immunity(self, rate = 0.87):
        '''immunity after t months have passed since exposure'''
        if self.t_infection > 30:   # !!! can t_infection be negative e.g. math.floor(-0.1) != 0
            self.current_immunity = max(1, self.postchallenge_peak_immunity * (self.t_infection/30) ** -rate) #t in this equation is measured in months, not days
        else:
            self.current_immunity = max(1, self.postchallenge_peak_immunity)
            
    def theta_Nab(self, a = 4.82,b =-0.30, c = 3.31, d = -0.32): 
        #prepopulated with the mle estimates from the mle model using the Gelfand data set- only biological theta is modeled
        Nab = self.prechallenge_immunity
        mean = a + b*np.log2(Nab) 
        stdev = np.sqrt(max(c + d*np.log2(Nab), 0))
        return np.exp(np.random.normal(loc=mean, scale=stdev))
        
    def calculate_naiive_shed_duration(self, params):
        w = sum([params.duration_s[locus] for locus in params.loci if (self.genome[locus]) & (locus not in ['syn', 'nonsyn', 'cpg', 'nonsyn_neutral'])])
        M = np.log(params.duration_params['u']) + w
        S = np.sqrt(np.absolute(1 + w))
            
        mu = M
        std = S * params.duration_params['sigma']
            

        duration = max(1, np.exp((mu) + np.random.normal(scale = std)))
        return duration
            
        
    
    def calculate_shed_duration(self,  mutation = None, params = config.params): #u= 12.5, delta = 1.16, sigma = 1.86
        """probability of shedding given Nab at time t (days post infection); 
        assumes that individual was infected at t = 0; time is measured in days
        Equation S1 in Famulaire paper
        delta_t = time (days) since last infection -- survival curve follows lognormal distribution"""
        
        if not mutation: #initialize
            w = sum([params.duration_s[locus] for locus in params.loci if (self.genome[locus]) & (locus not in ['syn', 'nonsyn', 'cpg', 'nonsyn_neutral'])])
            M = np.log(params.duration_params['u']) + w
            S = np.sqrt(np.absolute(1 + w))
            
            mu = M - params.duration_params['delta']*np.log2(self.prechallenge_immunity)
            std = S * params.duration_params['sigma']
            

            self.shed_duration = max(1, np.exp(mu + np.random.normal(scale = std)))
            
            self.naiive_shed_duration_reporter = max(1, np.exp(M + np.random.normal(scale = std)))
        else: #extend
            std = params.duration_params['sigma']
            self.shed_duration *= (np.exp(params.duration_s[mutation] + np.sqrt(np.absolute(params.duration_s[mutation])) * np.random.normal(scale = std)))
            self.shed_duration = max(1, self.shed_duration)
            self.naiive_shed_duration_reporter *= (np.exp(params.duration_s[mutation] + np.sqrt(np.absolute(params.duration_s[mutation])) * np.random.normal(scale = std)))
            self.naiive_shed_duration_reporter = max(1, self.naiive_shed_duration_reporter)
        
    def viral_shed(self, age, eta=1.65, v=0.17, epsilon =0.32):
        '''virus shed per gram, time is in months!'''
        if self.is_shed == True:
            predicted_concentration = 10**self.peak_cid50(age*12) * np.exp(eta - 0.5*v**2 - ((np.log(self.t_infection) - eta)**2) / (2*(v + epsilon*np.log(self.t_infection))**2) )/self.t_infection
            return max(10**2.6, predicted_concentration) 
        else: 
            return 0
        
        
    def peak_cid50(self, age, k = 0.056, Smax = 6.7, Smin = 4.3, tau = 12):
        '''returns the peak log10(cid50/g) given prior immunity, age is in months!'''
        if age >= 6:
            peak_cid50_naiive = (Smax - Smin)* np.exp((7 - age) / tau) + Smin
        else:
            peak_cid50_naiive = Smax
            
        return (1-k*np.log2(self.prechallenge_immunity))*peak_cid50_naiive
    
    
    def update(self, dt_years, params = config.params):
        #dt_years because the household model operates in terms of years
        if self.naiive != True:
            self.t_infection += dt_years * 365 #days!
            self.calculate_current_immunity()
            if self.t_infection >= self.shed_duration:
                self.is_shed = False

            if self.is_shed:
                for allele in self.genome:
                    if allele not in ['nonsyn', 'syn', 'cpg', 'nonsyn_neutral']:
                        if self.genome[allele] == 0:
                            if self.t_infection > self.mutation_wait_time[allele]:
                                self.genome[allele] = 1
                                self.calculate_shed_duration(allele)
                    else:
                    #elif allele != 'nonsyn_neutral':
                        if self.t_infection >= self.mutation_wait_time[allele]:
                            self.genome[allele] += 1
                            self.mutation_wait_time[allele] += self.mutate(params.fixation_rates[allele])
                    #else:
                    #    if self.optimum:
                    #        if self.t_infection >= self.mutation_wait_time[allele]:
                    #            self.genome[allele]+=1
                    #            self.mutation_wait_time[allele] += self.mutate(params.fixation_rates[allele])
                        
   
                if self.t_infection >= self.mutation_wait_time['recomb']:
                    if self.genome['nonsyn'] > 0:
                        N = np.random.randint(0,self.genome['nonsyn'])
                        self.genome['nonsyn'] += -1 * N
                    self.mutation_wait_time['recomb'] += self.mutate(params.recomb_rate)
                
                #if not self.optimum:
                #    self.check_epistasis(params)
    
    
    @classmethod
    def initialize_noninfant_immunity(cls, age, trial, strain_type='S2', params = config.params):
        """assumes that the mOPV2 challenge data is reflective of all strain types"""
        n_months = int(age) * 12
        I= ImmunoInfection(strain_type)
        current_immunity, current_ages = [], []                
        t_wait_time = np.random.gamma(shape = ImmunoInfection.infection_rate(0), scale =1)
        for i, t_step in enumerate(range(0, n_months)):
            current_age = 0 + (i+1) * 30/365
            I.update(30/365) #advance thirty days
            t_wait_time -= 1
            if t_wait_time < 0:
                I.attempt_infection(params = params)
                t_wait_time = np.random.gamma(shape = ImmunoInfection.infection_rate(current_age), scale =1)
        I.shed_duration = 0 #only use this to simulate immunity dynamics...? there is a low level of ongoing shedding using this method 
        I.is_shed = False
        I.genome = params.return_genome_template()
        return I
    
    @classmethod
    def initialize_S3_infant_immunity(cls, params = config.params):
        I = ImmunoInfection('S3')
        for day in range(18*7):
            I.update(1/365.)
            if day in 7*np.asarray([6,10,14]):
                I.attempt_infection(params=params) #posterior stdev = 0.032
        I.shed_duration = 0
        I.is_shed = False
        I.genome = params.return_genome_template()
        return I
    
    @classmethod
    def initialize_infant_immunity(cls, trial, params = config.params):
        '''S2 equivalent immunity following bOPV'''
        if trial == 'bOPV':
            I= ImmunoInfection('S1')
            for day in range(18*7):
                I.update(1/365.)
                if day in 7*np.asarray([6,10,14]):
                    I.attempt_infection(modifier = 0.3123757562833756, params = params) #posterior stdev = 0.032
            I.shed_duration = 0
            I.is_shed = False
            I.strain_type = 'S2'
            I.genome = params.return_genome_template()
            return I
        elif trial == 'tOPV':
            I= ImmunoInfection('S2')
            for day in range(18*7):
                I.update(1/365.)
                if day in 7*np.asarray([6,10,14]):
                    I.attempt_infection(modifier = 0.82, params = params) #posterior stdev = 0.0592141663766006
            I.shed_duration = 0
            I.is_shed = False
            I.genome = params.return_genome_template()
            return I
        else:
            raise ValueError('trial must be bOPV or tOPV')
    
    @classmethod
    def initialize_null_immunity(cls, strain_type = 'S2'):
        I = ImmunoInfection(strain_type)
        return I
        
    
    @classmethod
    def infection_rate(cls, age, alpha= 0.06085097476771715, beta = 24.624102979834937,  gamma = 2.2602043849994953):
        rate = beta * (1 - np.exp(-alpha*age) ) + gamma
        return rate
