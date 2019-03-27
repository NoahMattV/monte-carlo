# monte-carlo
ECE 845 Monte Carlo Scattering Simulation

Noah Van Der Weide

Bulk intrinsic GaAs. All electrons start in the gamma valley and are subject to
an electric field along the z-direction.

I make no claim that these results are entirely accurate and the main purpose of
this project is to showcase how the transport of particles in a Semiconductor
material can be modeled.

numOfParticles = number of particles under test
numOfTimeSteps = number of steps of time sized at deltaT
deltaT = time step size
Efield = electric field applied in z-direction.

It may be confusing, but in monte_carlo_top, E is for the kinetic energy of the
particle, whereas in generateScatt it represents the Efield. (E is still Ek in
getScattering). I did this because I'm irrational.
