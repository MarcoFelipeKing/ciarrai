# import libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime


# define constants
V = 36.247 # volume of air

# the analytical solution to the ODE
def conc_function(t,E,V,k, lambda_):
    return E/V * (1 - np.exp(-lambda_*t)+k*np.exp(-lambda_*t)) # K is a constant

# call a C++ function to solve E/V * (1 - np.exp(-lambda_*t)+k*np.exp(-lambda_*t)) for each time point

data = pd.read_csv('data/vmn_mp_conc.txt', sep='\t',header=0)

# convert time values to datetime objects
data["time"] = data["time"].apply(lambda x: datetime.strptime(x, "%M:%S"))

# calculate fraction of an hour by dividing the total number of minutes by 60
data["fraction_of_hour"] = data["time"].apply(lambda x: (x.minute+x.second/60) / 60)

#plot the analytical solution with the experimental data

# read in E_post, k_post, lambda_post
E_post = np.loadtxt('E_post.txt')
k_post = np.loadtxt('k_post.txt')
lambda_post = np.loadtxt('lambda_post.txt')

t0 = data['fraction_of_hour'][0:15]
t1 = data['fraction_of_hour'][15:43]
t2 = data['fraction_of_hour'][43:136]
#find mean of first 15 values from data['conc_08']
C0=data['conc_08'][0:15].mean()
#repeat C0 times to match length of t0
C0 = np.repeat(C0,len(t0))
C1 = conc_function(t1,E_post[0],V,0,lambda_post[0])
C2 = conc_function(t2,0,V,k_post[0],lambda_post[0])
C = np.concatenate((C0,C1,C2))
plt.plot(data['time'], data['conc_08'], 'o', label='Experimental data')
plt.plot(data['time'], C, label='Analytical solution')
plt.legend()
plt.show()