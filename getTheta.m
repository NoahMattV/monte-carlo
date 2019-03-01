% update theta

function out = getTheta()
  rtheta = 2*pi*rand; % between 0 and 2pi
  out = arccos(1-2*rtheta);
end
