import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
from tqdm import tqdm

def model(param,t):
    """
    Model function for the piecewise solution.
    
    Parameters
    ----------
    param : tuple
        Tuple containing the values of the parameters k and lambda.
    t : array-like
        Array of time values.
        
    Returns
    -------
    C : array-like
        Array of concentration values at the given times.
    """
    E, k, lambda_ = param
    
    # Initialize empty array for concentration values
    C = np.zeros_like(t)
    V = 36.247
    # Loop over time values
    for i, ti in enumerate(t):
        # If time is less than 5 minutes, set concentration to k
        if ti < 5/60:
            C[i] = k#0.003851195
        # If time is between 5 and 20 minutes, set concentration according to ODE solution
        elif 5/60 <= ti < 13/60:
            C[i] = k*np.exp(-lambda_*(ti-5/60)) + (E/V)*(1-np.exp(-lambda_*(ti-5/60)))/lambda_
        # If time is greater than 20 minutes, set concentration to 0
        else:
            C[i] = k*np.exp(-lambda_*(ti-13/60)) 
    # Add noise to the simulated data
    sigma = 0.001
    C += np.random.normal(0, sigma, size=C.shape)

    return C

def model_for_plotting(param,t):
    """
    Model function for the piecewise solution.
    
    Parameters
    ----------
    param : tuple
        Tuple containing the values of the parameters k and lambda.
    t : array-like
        Array of time values.
        
    Returns
    -------
    C : array-like
        Array of concentration values at the given times.
    """
    E, k, lambda_ = param
    
    # Initialize empty array for concentration values
    C = np.zeros_like(t)
    V = 36.247
    # Loop over time values
    for i, ti in enumerate(t):
        # If time is less than 5 minutes, set concentration to k
        if ti < 5/60:
            C[i] = 0.003851195
        # If time is between 5 and 20 minutes, set concentration according to ODE solution
        elif 5/60 <= ti < 13/60:
            C[i] = k*np.exp(-lambda_*(ti-5/60)) + (E/V)*(1-np.exp(-lambda_*(ti-5/60)))/lambda_
        # If time is greater than 20 minutes, set concentration to 0
        else:
            C[i] = k*np.exp(-lambda_*(ti-13/60)) 
    # Add noise to the simulated data
    #sigma = 0.001
    #C += np.random.normal(0, sigma, size=C.shape)

    return C

def distance(sim_data, data):
    """
    Distance function between simulated and observed data.
    
    Parameters
    ----------
    sim_data : array-like
        Simulated data.
    data : array-like
        Observed data.
        
    Returns
    -------
    d : float
        Distance between simulated and observed data.
    """
    # Calculate L2 norm between simulated and observed data
    d = np.linalg.norm(sim_data - data)
    
    return d

def prior(n_samples):
    """
    Prior distribution for the model parameters.
    
    Parameters
    ----------
    n_samples : int
        Number of samples to generate.
        
    Returns
    -------
    params : array-like
        Samples from the prior distribution.
    """
    # Generate samples from the prior distribution
    E = np.random.uniform(0, 10, n_samples)
    k = np.random.uniform(0, 10, n_samples)
    lambda_ = np.random.uniform(0, 10, n_samples)
    
    # Combine samples into array
    params = np.array([E, k, lambda_]).T
    
    return params



def abc(data, n_samples, epsilon, prior, model, distance):
    """
    Approximate Bayesian Computation (ABC) algorithm.
    
    Parameters
    ----------
    data : array-like
        Observed data.
    n_samples : int
        Number of samples to generate.
    epsilon : float
        Tolerance threshold for the distance between simulated and observed data.
    prior : callable
        Prior distribution for the model parameters.
    model : callable
        Model that generates simulated data.
    distance : callable
        Distance function between simulated and observed data.
        
    Returns
    -------
    samples : array-like
        Samples from the posterior distribution.
    """
    # Progress bar
    progress_bar = tqdm(total=n_samples)

    # Initialize empty array for samples
    samples = []
    
    # Generate samples from the prior distribution
    params = prior(n_samples)
    
    # Loop over samples
    for param in params:
        # Generate simulated data
        sim_data = model(param, np.linspace(0, 45/60, 135))
        
        # Compute distance between simulated and observed data
        d = distance(sim_data, data)

        # If distance is below threshold, add sample to array
        if d < epsilon:
            samples.append(param)
        
        progress_bar.update(1) 

    # order samples by distance d
    samples = np.array(samples)
    samples = samples[np.argsort(d)]

    progress_bar.close()   
    # Convert samples to numpy array and return
    return np.array(samples)

# Load data
data = pd.read_csv('data/vmn_mp_conc.txt', sep='\t',header=0)

# convert time values to datetime objects
data["time"] = data["time"].apply(lambda x: datetime.strptime(x, "%M:%S"))

# calculate fraction of an hour by dividing the total number of minutes by 60
data["fraction_of_hour"] = data["time"].apply(lambda x: (x.minute+x.second/60) / 60)

# subtract 0.003851195 from the concentration values
data["conc_2"] = data["conc_2"] #- 0.003851195

# extract time and concentration values
t = data["fraction_of_hour"].values
C = data["conc_2"].values

"""# Plot data
plt.plot(t, C, 'o', label='Data')
plt.xlabel('Time (h)')
plt.ylabel('Concentration (mug/m^3)')
plt.legend()
plt.show()"""

# Run ABC algorithm
samples = abc(C, 100_000, 0.125, prior, model, distance)

# Plot posterior distribution on 3 subplots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
axes[0].hist(samples[:, 0], bins=10, density=True)
axes[0].set_xlabel('E')
axes[1].hist(samples[:, 1], bins=10, density=True)
axes[1].set_xlabel('k')
axes[2].hist(samples[:, 2], bins=10, density=True)
axes[2].set_xlabel('lambda')
plt.show()

# Plot posterior predictive distribution
plt.plot(np.linspace(0, 45/60, 135), C, 'o', label='Data')
# plot the top value of the posterior predictive distribution
plt.plot(np.linspace(0, 45/60, 135), model_for_plotting(samples[np.argmax(samples[:, 0]), :],np.linspace(0, 45/60, 135)), label='Posterior predictive')
plt.xlabel('Time (h)')
plt.ylabel('Concentration (mug/m^3)')
#anotate the top value of the posterior predictive distributions on the graph
#plt.annotate('E = {:.2f}, k = {:.2f}, lambda = {:.2f}'.format(samples[np.argmax(samples[:, 0]), 0], samples[np.argmax(samples[:, 0]), 1], samples[np.argmax(samples[:, 0]), 2]), xy=(0.5, 0.5), xytext=(0.5, 0.5))#
plt.legend()
plt.show()

# print the top value of the posterior predictive distribution
print('E = {:.2f}, k = {:.2f}, lambda = {:.2f}'.format(samples[np.argmax(samples[:, 0]), 0], samples[np.argmax(samples[:, 0]), 1], samples[np.argmax(samples[:, 0]), 2]))
