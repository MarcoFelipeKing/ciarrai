# import libraries
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime


# define constants
V = 36.247 # volume of air

# the analytical solution to the ODE
def conc_function(t,E,V,k, lambda_):
    return E/V * (1 - np.exp(-lambda_*t))+k*np.exp(-lambda_*t) # K is a constant

# call a C++ function to solve E/V * (1 - np.exp(-lambda_*t)+k*np.exp(-lambda_*t)) for each time point

data = pd.read_csv('data/vmn_mp_conc.txt', sep='\t',header=0)

# convert time values to datetime objects
data["time"] = data["time"].apply(lambda x: datetime.strptime(x, "%M:%S"))

# calculate fraction of an hour by dividing the total number of minutes by 60
data["fraction_of_hour"] = data["time"].apply(lambda x: (x.minute+x.second/60) / 60)

#plot the analytical solution with the experimental data

# read in E_post, k_post, lambda_post
E_post = np.loadtxt('data/E_post.txt')
k_post = np.loadtxt('data/k_post.txt')
lambda_post = np.loadtxt('data/lambda_post.txt')

# plot the posterior distribution in 3 subplots
plt.figure(figsize=(5,9))
plt.subplot(3,1,1)
plt.hist(E_post, bins=20)
plt.xlabel('E')
plt.ylabel('Frequency')
plt.subplot(3,1,2)
plt.hist(k_post, bins=20)
plt.xlabel('k')
plt.ylabel('Frequency')
plt.subplot(3,1,3)
plt.hist(lambda_post, bins=20)
plt.xlabel('lambda')
plt.ylabel('Frequency')
plt.show()

t0 = data['fraction_of_hour'][0:15] # no emission between 0 and 5 minutes
t1 = data['fraction_of_hour'][15:44] # emission between 5 and 14 minutes
t2 = data['fraction_of_hour'][55:135]# no emission between 14 and 45 minutes
#find mean of first 15 values from data['conc_08']
C0=data['conc_2'][0:15].mean()
#repeat C0 times to match length of t0
C0 = np.repeat(C0,len(t0))
C1 = conc_function(t1,E_post[0],V,k_post[0],lambda_post[0])
C2 = conc_function(t2,0,V,k_post[0],lambda_post[0])
C = np.concatenate((C0,C1,C2))
_temp = np.array(data.loc[list(range(0,44)) + list(range(55, 135)), ['conc_2']])
_temp_time = np.array(data.loc[list(range(0,44)) + list(range(55, 135)), ['time']])
plt.plot(data['time'], data['conc_08'], 'o', label='80cm Experimental data')
plt.plot(data['time'], data['conc_2'], 'o', label='2m Experimental data')
plt.plot(_temp_time, C, label='Analytical solution')
#rotate x-axis labels
plt.xticks(rotation=15)
plt.xlabel('Time (min)')
plt.ylabel('Concentration (mg/m^3)')
plt.legend()
plt.show()



import pandas as pd 
df = pd.DataFrame({"E_post:":E_post,"k_post:":k_post,"lambda_post:":lambda_post})
print(df.describe())