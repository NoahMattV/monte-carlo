% Generate free-flight time for particle
function out = getTff();
  rt = rand;
  Gamma_o = 1; % I don't know what this is, so I'm setting it to 1.
  out = (-1/Gamma_o)*log(rt);
end
