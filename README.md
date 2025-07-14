Hybrid Monte-Carlo simulation code of (1+2)-d orbifold lattice Lagrangian used in arXiv:2506.00755. The Lagrangian is given by eq.(12) and eq.(13) in arXiv:2506.00755.

Fortran90 + Lapack. Compile main.f90.

You can specify the size of the gauge group, number of scalars, and number of lattice sites in size.h. Simulation parameters can be specified in input_v1.dat. 

Most of the input parameters should be self-explanatory, except for: 
- nskip (measurement is performed only when the number of trajectories in the HMC simulation is a multiple of nskip)
- ntau (number of steps N_{\tau} in HMC; see arXiv:1808.08490) and dtau_t and dtau_s (\Delta\tau for temporal and spatial directions (umat and zmat) in HMC; again, see arXiv:1808.08490)
- mersenne_seed (seed for the Mersenne Twister (random number generator)). Note that mersenne_seed is used only when init=1. 
