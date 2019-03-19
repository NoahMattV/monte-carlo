% update the kinetic energy of the particle

function [energy, energyInt] = getEnergy()

  energyInt = randi([1,1001], 1);
  %e = 1.602e-19; % already defined as a global variable, but kept in here just in case that doesn't work.
  E = linspace(1e-255,2*e,1001); % 0 to 2*e Joules (very close to 0).
  %re = 2*e*rand; % choose a number between 0 and 2*e
  out = -(3/2)*kT*log(re); % what to do with the -(3/2)*kT for the energyInt?
end
