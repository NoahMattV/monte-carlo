% update the kinetic energy of the particle

function out = getEnergy()
  e = 1.602e-19;
  re = 2*e*rand;
  out = -(3/2)*kT*log(re);
end
