% update the kinetic energy of the particle

function [energy, energyInt] = getEnergy()
  %e = 1.602e-19; % already defined as a global variable, but kept in here just in case that doesn't work.
  re = 2*e*rand; % choose a number between 0 and 2*e
  energy = -(3/2)*kT*log(re); % what to do with the -(3/2)*kT for the energyInt?
  energyInt = ceil(1001*energy); % integer value for energy to relate to arrays of scattering mechs.
end
