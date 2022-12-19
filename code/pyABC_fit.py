# Using pyabc to fit a model to data

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
from tqdm import tqdm
import pyabc as abc

import os
import tempfile

import matplotlib.pyplot as plt

from pyabc import ABCSMC, RV, Distribution, LocalTransition, MedianEpsilon
from pyabc.visualization import plot_data_callback, plot_kde_2d

db_path = "sqlite:///" + os.path.join(tempfile.gettempdir(), "test.db")

# measurement data
data = pd.read_csv('data/vmn_mp_conc.txt', sep='\t',header=0)
# convert time values to datetime objects
data["time"] = data["time"].apply(lambda x: datetime.strptime(x, "%M:%S"))

# calculate fraction of an hour by dividing the total number of minutes by 60
data["fraction_of_hour"] = data["time"].apply(lambda x: (x.minute+x.second/60) / 60)

data= np.array(data.loc[list(range(0,44)) + list(range(55, 135)), ['conc_2']])

# define constants
V = 36.247 # volume of air

def conc_function(t,E,V,k, lambda_):
    res = E/V * (1 - np.exp(-lambda_*t)) + k * np.exp(-lambda_*t)
    return  res

def deterministic_run(parameters):#precision,initial_contamination,r,C,d,g):
    #mp type
    t0 = data['fraction_of_hour'][0:15] # no emission between 0 and 5 minutes
    t1 = data['fraction_of_hour'][15:44] # emission between 5 and 14 minutes
    t2 = data['fraction_of_hour'][55:135]# no emission between 14 and 45 minutes
    C0 = data['conc_2'][0:15].mean() # find the mean of the first 15 data points
    C0 = np.repeat(C0,len(t0))
    V = 36.247
    C1 = conc_function(t1,parameters['E_prior'],V,parameters["k_prior"],parameters['lambda_prior']) # this is the mp experiment
    C2 = conc_function(t2,0,V,parameters["k_prior"],parameters['lambda_prior'])

    C = np.concatenate((C0,C1,C2))
    return{"Contamination":C }


def Distance(simulation, data):
    return np.absolute((data - simulation["Contamination"])).sum()

parameter_prior = Distribution(E_prior=RV("uniform", 0.0, 40.0),
                               k_prior=RV("uniform", 0.0, 10.0),
                               lambda_prior=RV("uniform", 0.0, 10.0))

#parameter_prior.get_parameter_names()

#No Noise
abc = ABCSMC(models=deterministic_run,
              parameter_priors=parameter_prior,
              distance_function=Distance,
              population_size=500,
              transitions=LocalTransition(k_fraction=.5),
              eps=MedianEpsilon(500, median_multiplier=0.7))

abc.new(db_path, {"Contamination": data})
history = abc.run(minimum_epsilon=2, max_nr_populations=10)