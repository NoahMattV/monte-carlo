% given the current valley, both the effective mass and scattering mechanism are
% determined. The effective mass of the electron is dependent on the final valley

function [s, v, m] = getScattering(v_in, Ek_int)
  % v_in = initial valley
  % Ek = kinetic energy of electron
  % s = scattering mechanism
  % v = final valley
  % m = effective mass after scattering

  % Scattering mechanism is chosen based on the energy of the electron.
  % The largest energy packet the electron can afford is the scattering mechanism employed.
  % The closest mechanism just above the electron energy (round up).

  % The new electron energy will most likely not be exactly equal to a discrete
  % energy in one of the scattering mechanisms, so rounding will have to occur for energy.
  % Professionals may look to see which step the energy is closest to and round up or down
  % accordingly, but I'm not getting paid to do this. Because I know how many steps there
  % are in the scattering mechanisms gammas, the getEnergy() function will assign an
  % integer in the range of that. To refer to the energy as it's value in Joules, I'll
  % have two variables for energy -- energy and energyInt. EnergyInt is the value that
  % relates to the scattering mech, whereas energy is the value converted to Joules,
  % according to the range used in scattering mechanism calculations
  % Between 0 and 2*e Joules with 1001 steps.

  % acoustic, pop absorption, pop emission, intervalley abs, intervalley em
  % from Gamma -> L (final valleys = 4), X (final valleys = 3)
  % from L -> L (final valleys = 3), X (final valleys = 3), Gamma (final valleys = 1)
  % from X -> X (final valleys = 2), L (final valleys = 4), Gamma (final vallyes = 1)
  % Final valleys based on degeneracy. L has 4 points of degeneracy, X has 3, Gamma has 1.
  global G_Tot = zeros(1001);

  v = randi([1,3], 1); % generates a random number in the range [1,3] inclusive (1, 2, or 3).
  % modify this to make sense in whatever respect you're using it in.

  if (v_in == 0) % Gamma
    s = randi([0,8], 1);
    switch s
      case 0
        % 0 = Acoustic
        G_Tot = G;
      case 1
        % 1 = Ac + POP abs
        G_Tot = G + Gpop_abs;

      case 2
        % 2 = Ac + POP abs + POP em
        G_Tot = G + Gpop_abs + Gpop_em;

      case 3
        % 3 = Ac + POP abs + POP em + Gamma to L abs (L = final valley)
        G_Tot = G + Gpop_abs + Gpop_em + Giv_abs_GL;
      case 4
        % 4 = Ac + POP abs + POP em + Gamma to L abs + Gamma to L em (L = final valley)
        G_Tot = G + Gpop_abs + Gpop_em + Giv_abs_GL + Giv_em_GL;

      case 5
        % 5 = Ac + POP abs + POP em + Gamma to X abs (X = final valley)
        G_Tot = G + Gpop_abs + Gpop_em + Giv_abs_GX;

      case 6
        % 6 = Ac + POP abs + POP em + Gamma to X abs + Gamma to X em (X = final valley)
        G_Tot = G + Gpop_abs + Gpop_em + Giv_abs_GX + Giv_em_GX;

      case 7
        % 7 = Self Scattering! (Do Nothing)

      otherwise
        % Self scattering
    end

  else if (v_in == 1) % X
    s = randi([0,10], 1);
    switch s
      case 0
        % 0 = Acoustic
        G_Tot = G;

      case 1
        % 1 = Ac + POP abs
        G_Tot = G + Gpop_abs;

      case 2
        % 2 = Ac + POP abs + POP em
        G_Tot = G + Gpop_abs + Gpop_em;

      case 3
        % 3 = Ac + POP abs + POP em + X to Gamma abs (Gamma = final valley)
        G_Tot = G + Gpop_abs + Gpop_em + Giv_abs_XG;

      case 4
        % 4 = Ac + POP abs + POP em + X to Gamma abs + X to Gamma em (Gamma = final valley)
        G_Tot = G + Gpop_abs + Gpop_em + Giv_abs_XG + Giv_em_XG;

      case 5
        % 5 = Ac + POP abs + POP em + X to L abs (L = final valley)
        G_Tot = G + Gpop_abs + Gpop_em + Giv_abs_XL;

      case 6
        % 6 = Ac + POP abs + POP em + X to L abs + X to L em (L = final valley)
        G_Tot = G + Gpop_abs + Gpop_em + Giv_abs_XL + Giv_em_XL;

      case 7
        % 7 = Ac + POP abs + POP em + X to X abs (X = final valley)
        G_Tot = G + Gpop_abs + Gpop_em + Giv_abs_XX;

      case 8
        % 8 = Ac + POP abs + POP em + X to X abs + X to X em (X = final valley)
        G_Tot = G + Gpop_abs + Gpop_em + Giv_abs_XX + Giv_em_XX;

      case 9
        % 9 = Self Scattering! (Do Nothing)

      otherwise
        % Self Scattering
    end

  else if (v_in == 2) % L
    s = randi([0,10], 1);
    s = randi([0,10], 1);
    switch s
      case 0
        % 0 = Acoustic
        G_Tot = G;
      case 1
        % 1 = Ac + POP abs
        G_Tot = G + Gpop_abs;

      case 2
        % 2 = Ac + POP abs + POP em
        G_Tot = G + Gpop_abs + Gpop_em;

      case 3
        % 3 = Ac + POP abs + POP em + L to Gamma abs (Gamma = final valley)
        G_Tot = G + Gpop_abs + Gpop_em + Giv_abs_LG;

      case 4
        % 4 = Ac + POP abs + POP em + L to Gamma abs + X to Gamma em (Gamma = final valley)
        G_Tot = G + Gpop_abs + Gpop_em + Giv_abs_LG + Giv_em_LG;

      case 5
        % 5 = Ac + POP abs + POP em + L to L abs (L = final valley)
        G_Tot = G + Gpop_abs + Gpop_em + Giv_abs_LL;

      case 6
        % 6 = Ac + POP abs + POP em + L to L abs + X to L em (L = final valley)
        G_Tot = G + Gpop_abs + Gpop_em + Giv_abs_LL + Giv_em_LL;

      case 7
        % 7 = Ac + POP abs + POP em + L to X abs (X = final valley)
        G_Tot = G + Gpop_abs + Gpop_em + Giv_abs_LX;

      case 8
        % 8 = Ac + POP abs + POP em + L to X abs + X to X em (X = final valley)
        G_Tot = G + Gpop_abs + Gpop_em + Giv_abs_LX + Giv_em_LX;

      case 9
        % 9 = Self Scattering! (Do Nothing)

      otherwise
        % Self Scattering
    end

  end

end
