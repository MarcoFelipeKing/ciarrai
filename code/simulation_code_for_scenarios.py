import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd
from scipy.integrate import odeint
from datetime import datetime

# make a measurement_times vector of every 20s for 1h and then convert to fraction of hour
measurement_times = np.arange(0, 3600, 20) / 3600

# Define model function
def model(param):
    """
    Model function for the piecewise solution.
    
    Parameters
    ----------
    param : dict
        Dictionary containing the values of the parameters E, lambda_1, and lambda_2.
    
    Returns
    -------
    C : dict
        Dictionary containing an array-like of concentration values at the corresponding measurement times.
    """
    
    E = param['E']
    lambda_1 = param['lambda_1']
    lambda_2 = param['lambda_2']
    V = param['V']
   
    t = measurement_times
    C0 = 0 #starting concentration of the data from first 15 readings (first 5 minutes)
    # Initialize empty array for concentration values
    C = np.zeros_like(t)
    
    
    # Loop over time values
    for i, ti in enumerate(t):
        # If time is less than 5 minutes, set concentration to constant background value
        if ti < 5/60:
            C[i] = C0
        # If time is between 5 and 14 minutes, set concentration according to ODE solution
        elif 5/60 <= ti < 14/60:
            C_at_5_minutes = C0
            C[i] = (E/V)*(1-np.exp(-lambda_1*(ti-5/60)))/lambda_1 + C_at_5_minutes * np.exp(-lambda_1*(ti-5/60))
            C_14=C[i] #last concentration value at 14 minutes
        # If time is between 14 and 17 minutes, set concentration according to fast decay part
        elif 14/60 <= ti < 17/60:
            C[i] = C_14 * np.exp(-lambda_1*(ti-14/60))
            C_17=C[i] #last concentration value at 17 minutes
        # If time is greater than 17 minutes, set concentration according to slow decay part
        else:
            C[i] = C_17 * np.exp(-lambda_2*(ti-17/60))

    # Add noise to the simulated data
    sigma = 0.001
    C += np.random.normal(0, sigma, size=C.shape)

    return {"Concentration": C}

import numpy as np
from scipy.integrate import cumtrapz

# Define function to calculate dose
def calculate_dose(C, p, dt):
    """
    Calculate the total dose inhaled over time.
    Parameters:
    C: numpy array of drug concentration in the air over time
    p: breathing volume/hour
    dt: time step size
    Returns:
    D: numpy array of total dose inhaled over time
    """
    dD_dt = p * C
    D = cumtrapz(dD_dt, dx=dt, initial=0)
    return D

# Calculate inhaled dose for healthcare worker at 2.0 m for 10 minutes
#C = model(param)["Concentration"]
# p = 0.54, breathing rate (m3/hr) according to EPA exposure factor guidelines for Sedentary, defined as sitting and standing (see Table 6-40 for original data).
# p = 1.38, breathing rate (m3/hr) according to EPA exposure factor guidelines for Light activity, defined as walking at speed level 1.5−3.0 mph (see Table 6-40 for original data).

p = 1.38
dt = 20 / 3600
#new_C = C[14:44]
#dose_ = calculate_dose(new_C,p,dt)
#take last value of dose_ to get total dose inhaled over time
#dose_ = dose_[-1]

"""
*changes made*

Across all nebulisers and interfaces tested, on average, 
lambda_1 = 6.5 hr-1 and 
lambda_2 = 1.3 hr-1

Assuming the AER increases by 0.55 hr-1 for opening one window
and increases additionally by 0.55 hr-1 when a door is opened,

The following lambda values apply

Scenario_1 : No window or door open
            lambda_1 = 6.5 hr-1 
            lambda_2 = 1.3 hr-1
            
            Workings:
            lambda_1 = 6.5 hr-1 
            (AER + Deposition Rate (DR) + Mixing factor) 
            (0.55 + 5.95 + 0) #mixing factor goes to zero 
            
            lambda_2 = 1.3 hr-1
            (0.55 + 0.75 + 0)  
            
Scenario_2 : Window open
            lambda_1 = 7.05 hr-1
            lambda_2 = 1.85 hr-1
            
            (add 0.55 hr-1 to account for open window)
            
Scenario_2 : Window and door open
            lambda_1 = 7.6 hr-1
            lambda_2 = 2.4 hr-1
            
            (add 0.55 hr-1 to account for open door)            
            
How much does DR decrease when AER is increased?
"""

#Create an array of scenario_1 and scenario_2 parameters to feed into the model


# Define scenarios
scenarios = {
    "JN_ad_mp": {"E": 141.33, "lambda_1": 10.56, "lambda_2": 1.22, "V": 28.00},
    "JN_ad_fm": {"E": 259.73, "lambda_1": 13.58, "lambda_2": 1.04, "V": 28.00},
    "JN_pd_fm": {"E": 72.03, "lambda_1": 1.82, "lambda_2": 1.01, "V": 28.00},
    "VMN_ad_mp": {"E": 10.95, "lambda_1": 8.59, "lambda_2": 1.64, "V": 28.00},
    "VMN_ad_fm": {"E": 18.90, "lambda_1": 3.53, "lambda_2": 1.49, "V": 28.00},
    "VMN_pd_mp": {"E": 10.11, "lambda_1": 6.35, "lambda_2": 1.74, "V": 28.00},
    "VMN_pd_fm": {"E": 12.34, "lambda_1": 1.04, "lambda_2": 1.01, "V": 28.00},
}

# Initialize empty DataFrame to store results
results = pd.DataFrame(columns=["Scenario", "Total Dose"])

# Loop over scenarios
for scenario, param in scenarios.items():
    # Calculate concentration
    C = model(param)["Concentration"]
    # Calculate dose
    new_C = C[14:44] #only take values between 14 and 44 minutes, is that what you want?
    dose_ = calculate_dose(new_C, p, dt)
    total_dose = dose_[-1]
    # Append results to DataFrame
    results = results.append({"Scenario": scenario, "Total Dose": total_dose}, ignore_index=True)

# Plot results
plt.figure(figsize=(10, 6))
plt.bar(results["Scenario"], results["Total Dose"])
plt.xlabel("Scenario")
plt.ylabel("Total Dose")
plt.title("Total Dose for Different Scenarios")
plt.xticks(rotation=45)
plt.show()