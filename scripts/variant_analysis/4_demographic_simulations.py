#!/usr/bin/env python

#Run a population simulation

#NB written for Python 3.5+

################
# DEPENDENCIES #
################

import os
import pickle
import pandas as pd
import numpy as np
import simuPOP as sim

import serial_functions as serf #custom functions (./serial_functions.py)

#Define path to interim data
INTERIM_PATH = os.path.abspath(os.path.join(CWD, os.pardir, os.pardir, 'data', 'interim'))

#set up the simulation
#NOTE THAT THIS OUTCOME WILL BE DIFFERENT EVERY TIME YOU RUN THE SIMULATION

simulation_outcome = {}

#Define sampling rate:
for sampling in [36,48,60,72,84]:

    #Define the population
    pop = sim.Population(size=100, loci=1)

    #Define the demographic trajectory
    demo_function = serf.demo_dynamic(101, 100., np.log(2), 100000., 10., 30, sampling)

    #Initiate the simulation
    simu = sim.Simulator(pop,rep=100)

    #Evolve population
    simu.evolve(
        initOps=[
            sim.InitSex(),
            sim.InitGenotype(freq=[0.95,0.05]), #proportion of minor to major allele
            sim.PyExec('traj=[]') #Initiate a recipient list for frequency outputs
        ],
        matingScheme=sim.RandomSelection(subPopSize=demo_function), #random binary fission
        postOps=[
            sim.Stat(alleleFreq=0),
            sim.PyExec('traj.append(alleleFreq[0][1])', step=10), #record state of the simulation
        ],
        gen=101 #number of generations over which to run the simulation
    )

    #Store simulation output, specifically the state of the alleles at the
    #last recorded simulated timepoint.
    simulation_outcome[sampling]=[simu.dvars(x).traj[-1] for x in range(100)]

#Reformat the simulation outcome
simulation_outcome = pd.DataFrame(simulation_outcome)

#Save outcome to csv.
simulation_outcome.to_csv('{}/6_simuPOP_outcome.csv'.format(INTERIM_PATH))
