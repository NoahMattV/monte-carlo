% Generate free-flight time for particle
function out = getTff()
  rt = rand;

  global GTot_X_Max;
  Gamma_o = max(GTot_X_Max);
  out = (-1/Gamma_o)*log(rt);
  %out = 50;
end
