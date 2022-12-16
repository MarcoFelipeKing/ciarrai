import numpy as np

def point_source_pde(C, D, k, dx, dt):
    """Solves the PDE for a continuous point source under diffusion and decay.
    
    Parameters:
    C (ndarray): Concentration at each location x at time t.
    D (float): Diffusion coefficient.
    k (float): Decay rate.
    dx (float): Spatial discretization step.
    dt (float): Time discretization step.
    
    Returns:
    ndarray: Concentration at each location x at time t + dt.
    """
    # Number of grid points
    N = C.shape[0]
    
    # Initialize array for updated concentration values
    C_new = np.zeros(N)
    
    # Iterate over all grid points
    for i in range(N):
        # Implement finite difference formula for diffusion
        C_diffusion = D * (C[(i+1)%N] - 2*C[i] + C[(i-1)%N]) / dx**2
        
        # Implement finite difference formula for decay
        C_decay = -k * C[i]
        
        # Update concentration at each grid point
        C_new[i] = C[i] + dt * (C_diffusion + C_decay)
    
    return C_new


# Set diffusion coefficient and decay rate
D = 0.1
k = 0.01

# Set spatial and temporal discretization steps
dx = 0.1
dt = 0.001

# Set initial concentration
C = np.zeros(100)
C[50] = 1

# Solve PDE for 1000 time steps
for t in range(1000):
    C = point_source_pde(C, D, k, dx, dt)
#print(C)
# Plot concentration profile at 5 different times
import matplotlib.pyplot as plt
plt.plot(C)
plt.xlabel('Distance [x]')
plt.ylabel('C')
plt.show()
