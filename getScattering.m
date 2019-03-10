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


  v = randi([1,3], 1); % generates a random number in the range [1,3] inclusive (1, 2, or 3).
  % modify this to make sense in whatever respect you're using it in.

  if (v_in == 0) % Gamma
    s = randi([0,8], 1);
    % 0 = Acoustic

    % 1 = Ac + II
    % 2 = Ac + II + POP abs
    % 3 = Ac + II + POP abs + POP em
    % 4 = Ac + II + POP abs + POP em + Gamma to L abs (L = final valley)
    % 5 = Ac + II + POP abs + POP em + Gamma to L abs + Gamma to L em (L = final valley)
    % 6 = Ac + II + POP abs + POP em + Gamma to X abs (X = final valley)
    % 7 = Ac + II + POP abs + POP em + Gamma to X abs + Gamma to X em (X = final valley)
    % 8 = Self Scattering! (Do Nothing)

  else if (v_in == 1) % X
    s = randi([0,10], 1);
    switch s
      case 0
      % 0 = Acoustic

      case 1
      % 1 = Ac + II
      % 2 = Ac + II + POP abs
      % 3 = Ac + II + POP abs + POP em
      % 4 = Ac + II + POP abs + POP em + X to Gamma abs (Gamma = final valley)
      % 5 = Ac + II + POP abs + POP em + X to Gamma abs + X to Gamma em (Gamma = final valley)
      % 6 = Ac + II + POP abs + POP em + X to L abs (L = final valley)
      % 7 = Ac + II + POP abs + POP em + X to L abs + X to L em (L = final valley)
      % 8 = Ac + II + POP abs + POP em + X to X abs (X = final valley)
      % 9 = Ac + II + POP abs + POP em + X to X abs + X to X em (X = final valley)
      % 10 = Self Scattering! (Do Nothing)

  else if (v_in == 2) % L
    s = randi([0,10], 1);

    % 0 = Acoustic
    % 1 = Ac + II
    % 2 = Ac + II + POP abs
    % 3 = Ac + II + POP abs + POP em
    % 4 = Ac + II + POP abs + POP em + L to Gamma abs (Gamma = final valley)
    % 5 = Ac + II + POP abs + POP em + L to Gamma abs + L to Gamma em (Gamma = final valley)
    % 6 = Ac + II + POP abs + POP em + L to L abs (L = final valley)
    % 7 = Ac + II + POP abs + POP em + L to L abs + L to L em (L = final valley)
    % 8 = Ac + II + POP abs + POP em + L to X abs (X = final valley)
    % 9 = Ac + II + POP abs + POP em + L to X abs + X to X em (X = final valley)
    % 10 = Self Scattering! (Do Nothing)
  end

end
