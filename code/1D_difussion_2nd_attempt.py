""" The 1D diffusion equation with a decay term can be written as:

$$\frac{\partial C}{\partial t} = D \frac{\partial^2 C}{\partial x^2} - kC$$

where C is the concentration of the contaminant, t is time, x is position, D is the diffusion coefficient, and k is the decay rate.

To solve this equation, we will need to specify the initial and boundary conditions.

Initial condition: At time t=0, the concentration of the contaminant is C(x,0) = C0(x).

Boundary conditions: Since there are no doors or windows in the room, we can assume that the contaminant is confined to the room. Therefore, the concentration of the contaminant at the walls of the room (x=0 and x=L) is zero: C(0,t) = C(L,t) = 0.

To solve this equation, we can use the finite difference method, which involves discretizing the space and time dimensions and approximating the derivatives with finite differences.

We can start by discretizing the space dimension into N equally spaced points, with a spacing of dx = L/(N-1). The time dimension can be discretized into M equally spaced points, with a spacing of dt = T/M, where T is the total time for which the contaminant is released.

The finite difference approximation of the second derivative in the x-direction is:

$$\frac{\partial^2 C}{\partial x^2} \approx \frac{C_{i+1} - 2C_i + C_{i-1}}{dx^2}$$

where C_i is the concentration at position x_i = i*dx.

Substituting this into the diffusion equation, we get:

$$\frac{C_i^{m+1} - C_i^m}{dt} = D \frac{C_{i+1}^m - 2C_i^m + C_{i-1}^m}{dx^2} - kC_i^m$$

where m is the time index and the superscript refers to the value at time t_m = m*dt.

We can rearrange this equation to solve for C_i^{m+1}, the concentration at the next time step:

$$C_i^{m+1} = C_i^m + D \frac{dt}{dx^2} (C_{i+1}^m - 2C_i^m + C_{i-1}^m) - kdtC_i^m$$

This equation can be used to update the concentration at each grid point at each time step, using the values of the concentration at the previous time step.

To solve the equation, we can start by specifying the initial concentration C0(x) and the diffusion coefficient D, decay rate k, and the total time T for which the contaminant is released. We can then iterate over time steps, updating the concentration at each grid point using the finite difference equation above, until we reach the final time T.

I hope this helps! Let me know if you have any questions or need further clarification."""


import numpy as np

# Set the diffusion coefficient, decay rate, and total time
D = 1.0
k = 0.1
T = 100

# Set the number of grid points and the grid spacing
N = 1001
L = 1.0
dx = L/(N-1)

# Set the time step size
dt = 0.001
M = int(T/dt)

# Set the initial concentration
C0 = np.zeros(N)
C0[int(N/2)] = 1.0

# Set the boundary conditions
C0[0] = 0.0
C0[-1] = 0.0

# Set the initial concentration array
C = np.zeros((M+1,N))
C[0,:] = C0

# Iterate over time steps
for m in range(1,M+1):
    # Set the source term
    S = np.zeros(N)
    if m*dt <= 10.0:
        S[int(N/2)] = 1.0/dt
    # Update the concentration at each grid point
    for i in range(1,N-1):
        C[m,i] = C[m-1,i] + D*dt/dx**2*(C[m-1,i+1] - 2*C[m-1,i] + C[m-1,i-1]) - k*dt*C[m-1,i] + S[i]*dt
    # Enforce the boundary conditions
    C[m,0] = 0.0
    C[m,-1] = 0.0

# Plot the concentration profile at the final time
import matplotlib.pyplot as plt
plt.plot(C[-1,:])
plt.show()

$$C_i^{m+1} = \frac{C_i^m + C_i^{m+1} + D \frac{dt}{2dx^2} (C_{i+1}^m - 2C_i^m + C_{i-1}^m) - kdtC_i^m}{1 + kdt}$$