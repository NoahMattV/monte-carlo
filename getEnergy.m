% update the kinetic energy of the particle

function [energy, energyInt] = getEnergy()
  global e;
  global kT;
  E = linspace(1e-255,2*e,1001); % 0 to 2*e Joules
  re = 2*e*rand; % choose a number between 0 and 2*e
  %re = rand;
  energy = -(3/2)*kT*log(re); % what to do with the -(3/2)*kT for the energyInt?
  %energyInt = abs(floor(1001*(energy/(2*e)))); % integer value for energy to relate to arrays of scattering mechs.
  A = repmat(energy, [1 length(E)]);
  [minValue,energyInt] = min(abs(A-E));
  %closestValue = energy(closestIndex);

end
