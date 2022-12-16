# Solve the ODE for airborne particle conentration
# dC/dt = E/V - lambda*C
# C = concentration 
# E = emission rate
# V = volume of air
# lambda = decay rate (deposition and ventilation removal)

# We will compare the results of the ODE solution with experimental data to infer the decay rate and emission rate

# Looking at the data, we can see that the emission rate is constant between 5 and 14 minutes but the decay is not well described by a single decay rate
# We will use an ABC algorithm to infer the decay rate and emission rate


# import libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
from tqdm import tqdm

# define constants
V = 36.247 # volume of air

# the analytical solution to the ODE
def conc_function(t,E,V,k, lambda_):
    return E/V * (1 - np.exp(-lambda_*t)) + k*np.exp(-lambda_*t) # K is a constant

# call a C++ function to solve E/V * (1 - np.exp(-lambda_*t)+k*np.exp(-lambda_*t)) for each time point

data = pd.read_csv('data/vmn_mp_conc.txt', sep='\t',header=0)
data_valved = pd.read_csv('data/vmn_ad_conc.txt', sep='\t',header=0)

# convert time values to datetime objects
data["time"] = data["time"].apply(lambda x: datetime.strptime(x, "%M:%S"))
data_valved["time"] = data['time']#data_valved["time"].apply(lambda x: datetime.strptime(x, "%M:%S"))

# calculate fraction of an hour by dividing the total number of minutes by 60
data["fraction_of_hour"] = data["time"].apply(lambda x: (x.minute+x.second/60) / 60)
data_valved["fraction_of_hour"] = data['fraction_of_hour']#data_valved["time"].apply(lambda x: (x.minute+x.second/60) / 60)
# plot the experimental data
#plt.plot(data['time'], data['conc_08'], 'o', label='Experimental data')

# program an ABC algorithm to infer the decay rate lambda_, emission rate E and constant k

N=1000

# define the distance function
# the distance function is the sum of the squared differences between the model and the data
def distance(model,data):
    return np.sum((model - data)**2)

# define the ABC algorithm
# for 1000 different E and k values run the model and calculate the distance
# if the distance is less than 0.1, save the E and k values
# Initialize the progress bar
def ABC():
    print('Running ABC algorithm...')
    progress_bar = tqdm(total=N)
    E_post = []
    k_post = []
    lambda_post = []
    dist = []

    #mp type
    t0 = data['fraction_of_hour'][0:15] # no emission between 0 and 5 minutes
    t1 = data['fraction_of_hour'][15:44] # emission between 5 and 14 minutes
    t2 = data['fraction_of_hour'][55:135]# no emission between 14 and 45 minutes
    C0 = data['conc_2'][0:15].mean() # find the mean of the first 15 data points
    C0 = np.repeat(C0,len(t0)) # repeat the mean value to match the length of t0
    #C0_2m = np.repeat(data['conc_2'][0:15].mean(),len(t0))
    """
    #valved 
    t0_valved = data_valved['fraction_of_hour'][0:15] # no emission between 0 and 5 minutes
    t1_valved = data_valved['fraction_of_hour'][15:46] # emission between 5 and 14 minutes
    t2_valved = data_valved['fraction_of_hour'][46:136]# no emission between 14 and 45 minutes
    C0_valved = data_valved['conc_08'][0:15].mean() # find the mean of the first 15 data points
    C0_valved = np.repeat(C0_valved,len(t0_valved)) # repeat the mean value to match the length of t0
    C0_2m_valved = np.repeat(data_valved['conc_2'][0:15].mean(),len(t0_valved)) """

    i=0 # iteration counter
    #accepted = 0 # accepted counter
    epsilon = 20# distance threshold
    while len(E_post) < N:
        i+=1
        k_prior = np.random.uniform(1E-1,1,1)
        E_prior = np.random.uniform(1E-2,50,1)
        lambda_prior = np.random.uniform(1E-2,3,1)
        # solve piece-wise. between 0 and 10 minutes, E= prior value
        # between 0 and 5 minutes E = 0. between 5 and 14 minutes E = prior value. between 14 and 45 minutes E = 0
        # solve C at the same time points as the data['time']
        # convert data['time'] to time then to fraction of an hour
        #repeat C0 times to match length of t0
        #MP type 
        C1 = conc_function(t1,E_prior,V,0,29) # this is the mp experiment
        C2 = conc_function(t2,0,V,k_prior,lambda_prior)
        # add Gaussian noise to the mode
        C0 = C0 #+ np.random.normal(0,0.005,len(C0))
        C1 = C1 + np.random.normal(0,0.01,len(C1))
        C2 = C2 + np.random.normal(0,0.005,len(C2))
        C = np.concatenate((C0,C1,C2))

        """C1_2m = conc_function(t1,E_prior,V,0,2*lambda_prior)
        C2_2m = conc_function(t2,0,V,k_prior,2*lambda_prior)
        C_2m = np.concatenate((C0_2m,C1_2m,C2_2m))

        # valved type
        C1_valved = conc_function(t1_valved,E_prior,V,0,lambda_prior)
        C2_valved = conc_function(t2_valved,0,V,k_prior,lambda_prior)
        # add Gaussian noise to the mode
        C0_valved = C0_valved #+ np.random.normal(0,0.01,len(C0_valved))
        C1_valved = C1_valved #+ np.random.normal(0,0.01,len(C1_valved))
        C2_valved = C2_valved #+ np.random.normal(0,0.01,len(C2_valved))
        C_valved = np.concatenate((C0_valved,C1_valved,C2_valved))

        C1_2m_valved = conc_function(t1_valved,E_prior,V,0,2*lambda_prior)
        C2_2m_valved = conc_function(t2_valved,0,V,k_prior,2*lambda_prior)
        C_2m_valved = np.concatenate((C0_2m_valved,C1_2m_valved,C2_2m_valved))"""

#select 0:15, 15:43, 55:136
        _temp = np.array(data.loc[list(range(0,44)) + list(range(55, 135)), ['conc_2']])
        #print(_temp)
        #print(C)
        delta = distance(C,_temp)#+distance(C_2m,data['conc_2'])+distance(C_valved,data_valved['conc_08'])+distance(C_2m_valved,data_valved['conc_2'])
        if delta < epsilon:
            E_post.append(E_prior)
            k_post.append(k_prior)
            lambda_post.append(lambda_prior)
            dist.append(delta)
            #accepted+=1
            progress_bar.update(1)      
    progress_bar.close()
    #print the acceptance rate
    print('Final Acceptance:',N/i*100,'%')
    return E_post, k_post, lambda_post, dist


# run the ABC algorithm
E_post, k_post,lambda_post, dist = ABC()

# function to sort the posterior distributions by distance
def sort_posterior(E_post,k_post,lambda_post,dist):
    dist=np.array(dist)
    idx = dist.argsort()
    E_post = np.array(E_post)
    k_post = np.array(k_post)
    lambda_post = np.array(lambda_post)

    E_post = E_post[idx]
    k_post = k_post[idx]
    lambda_post = lambda_post[idx]
    dist = dist[idx]
    return E_post, k_post,lambda_post,dist

# sort the posterior distributions by distance
E_post, k_post,lambda_post,dist = sort_posterior(E_post,k_post,lambda_post,dist)

#save the posterior distributions
np.savetxt('data/E_post.txt', E_post)
np.savetxt('data/k_post.txt', k_post)
np.savetxt('data/lambda_post.txt', lambda_post)
np.savetxt('data/dist.txt', dist)
