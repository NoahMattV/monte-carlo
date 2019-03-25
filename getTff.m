% Generate free-flight time for particle
function out = getTff()
  rt = rand;

  global GTot_L_Max;
  Gamma_o = max(GTot_L_Max);
  out = (-1/Gamma_o)*log(rt);
end
