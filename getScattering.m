% given the current valley, both the effective mass and scattering mechanism are
% determined. The effective mass of the electron is dependent on the final valley

function [s, v, m] = getScattering(v_in, Ek)
  % v_in = initial valley
  % Ek = kinetic energy of electron
  % s = scattering mechanism
  % v = final valley
  % m = effective mass after scattering

  % acoustic, ionized impurity, pop absorption, pop emission, intervalley abs, intervalley em
  % from Gamma -> L (final valleys = 4), X (final valleys = 3)
  % from L -> L (final valleys = 3), X (final valleys = 3), Gamma (final valleys = 1)
  % from X -> X (final valleys = 2), L (final valleys = 4), Gamma (final vallyes = 1)
  % Final valleys based on degeneracy. L has 4 points of degeneracy, X has 3, Gamma has 1.

  

end
