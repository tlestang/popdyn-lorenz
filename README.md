# popdyn-lorenz
Population dynamics with Lorenz model simulated through the Euler scheme.
Parameter for both the population dynamics and Lorenz model must, for now, be directly inputed in the code. If you do, please recompile the code right after that by typing 'make' in the command line.

Parameters for the population dynamics can be found at the very beginning of int main() in popdyn.cpp.
T : Total length of the file serie.
K : Number of realisations.
Nc : Fixed number of clones.
dT : Time between two cloning events.

Parameters for the Lorenz model can be found at the very beginning of void lorenz_transient() and double lorenz() in their respective .cpp file. MAKE SURE THESE PARAMETERS ARE THE SAME IN BOTH FILES.
timeStep : timestep for Euler scheme
sigma, rho, beta : parameters for the Lorenz equations of movement.
double *x0 : 3-component vector containing initial conditions.
double *m_x : final 3-component vector containing the values for the coordinates at the end of the iterations. This is needed because intermediate coordinate are not stored.
