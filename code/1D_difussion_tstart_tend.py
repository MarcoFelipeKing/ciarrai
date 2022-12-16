# 1D diffusion with a piece-wise release of contaminant

import numpy as np
def point_source_pde(C, D, k, dx, dt, t, t_start, t_end):
    """Solves the PDE for a continuous point source under diffusion and decay, with a continuous release of contaminant from t_start to t_end.
    
    Parameters:
    C (ndarray): Concentration at each location x at time t.
    D (float): Diffusion coefficient.
    k (float): Decay rate.
    dx (float): Spatial discretization step.
    dt (float): Time discretization step.
    t (int): Current time step.
    t_start (int): Start time of contaminant release.
    t_end (int): End time of contaminant release.
    
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
        
        # Implement continuous release of contaminant at center
        if t_start <= t < t_end:
            C_release = 1
        else:
            C_release = 0
        
        # Update concentration at each grid point
        C_new[i] = C[i] + dt * (C_diffusion + C_decay + C_release)
    
    return C_new


# Set spatial and temporal discretization steps
dx = 0.01
dt = 0.001

# Set release start and end times
t_start = 5
t_end = 15

# Set initial concentration
C = np.zeros(100)

# Solve PDE for release period
for t in range(t_start, t_end):
    C = point_source_pde(C, D, k, dx, dt, t, t_start, t_end)

# Solve PDE for 1000 time steps after the release period
for t in range(t_end, t_end+1000):
    C = point_source_pde(C, D, k, dx, dt, t, t_start, t_end)
