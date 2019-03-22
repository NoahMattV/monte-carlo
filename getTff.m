% Generate free-flight time for particle
function out = getTff()
  rt = rand;
  %Gamma_o = 1; % I don't know what this is, so I'm setting it to 1.
  global GTot_L_Max;
  Gamma_o = max(GTot_L_Max);
  out = (-1/Gamma_o)*log(rt);
end
