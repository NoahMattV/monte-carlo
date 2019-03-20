% given the current valley, both the effective mass and scattering mechanism are
% determined. The effective mass of the electron is dependent on the final valley

function [s, v, m] = getScattering(v_in, E_int)
  % v_in = initial valley
  % Ek = kinetic energy of electron
  % s = scattering mechanism (1 = ac, 2 = pop abs, 3 = pop em, 4 = iv abs, 5 = iv em, 6 = self)
  % v = final valley (1 = Gamma, 2 = X, 3 = L)
  % m = effective mass after scattering

  % acoustic, pop absorption, pop emission, intervalley abs, intervalley em
  % from Gamma -> L (final valleys = 4), X (final valleys = 3)
  % from L -> L (final valleys = 3), X (final valleys = 3), Gamma (final valleys = 1)
  % from X -> X (final valleys = 2), L (final valleys = 4), Gamma (final vallyes = 1)
  % Final valleys based on degeneracy. L has 4 points of degeneracy, X has 3, Gamma has 1.

  global G_Tot = zeros(1001);
  m0 = 9.11e-31 %kg
  rg = 0;
  v = randi([1,3], 1); % generates a random number in the range [1,3] inclusive (1, 2, or 3).
  % modify this to make sense in whatever respect you're using it in.
  r = rand();

  switch v_in
    case 1 % Gamma
      rg = r*GTot_G_Max;
      if (rg <= GTot_G(1, E_int))
        % acoustic (elastic, isotropic)
        s = 1;
        v = v_in;
        m = 0.067*m0;

      elseif (rg <= GTot_G(2, E_int))
        % POP abs (inelastic, anisotropic)
        s = 2;
        v = v_in;
        m = 0.067*m0;

      elseif (rg <= GTot_G(3, E_int))
        % POP em (inlastic, anisotropic)
        s = 3;
        v = v_in;
        m = 0.067*m0;

      elseif (rg <= GTot_G(4, E_int))
        % Intervalley G to L abs (inelastic, isotropic)
        s = 4;
        v = 3;
        m = 0.85*m0; % effective mass of density of states m = (16*ml*mt^2)^(1/3)

      elseif (rg <= GTot_G(5, E_int))
        % Intervalley G to L em (inelastic, isotropic)
        s = 5;
        v = 3;
        m = 0.85*m0;

      elseif (rg <= GTot_G(6, E_int))
        % Intervalley G to X abs (inelastic, isotropic)
        s = 4;
        v = 2;
        m = 0.85*m0; % effective mass of density of states m = (9*ml*mt^2)^(1/3)

      elseif (rg <= GTot_G(7, E_int))
        % Intervalley G to X em (inelastic, isotropic)
        s = 5;
        v = 2;
        m = 0.85*m0;

      else
        % Self-scattering -- do nothing
        s = 6;
        v = v_in;
        m = 0.067;
      end

    case 2 % X
      rg = r*GTot_X_Max;
      if (rg <= GTot_X(1, E_int))
        % acoustic (elastic, isotropic)
        s = 1;
        v = v_in;
        m = 0.85*m0;

      elseif (rg <= GTot_X(2, E_int))
        % POP abs (inelastic, anisotropic)
        s = 2;
        v = v_in;
        m = 0.85*m0;

      elseif (rg <= GTot_X(3, E_int))
        % POP em (inelastic, anisotropic)
        s = 3;
        v = v_in;
        m = 0.85*m0;

      elseif (rg <= GTot_X(4, E_int))
        % Intervalley X to G abs (inelastic, isotropic)
        s = 4;
        v = 1;
        m = 0.067*m0;
      elseif (rg <= GTot_X(5, E_int))
        % Intervalley X to G em (inelastic, isotropic)
        s = 5;
        v = 1;
        m = 0.067*m0;

      elseif (rg <= GTot_X(6, E_int))
        % Intervalley X to L abs (inelastic, isotropic)
        s = 4;
        v = 3;
        m = 0.85*m0;

      elseif (rg <= GTot_X(7, E_int))
        % Intervalley X to L em (inelastic, isotropic)
        s = 5;
        v = 3;
        m = 0.85*m0;

      elseif (rg <= GTot_X(8, E_int))
        % Intervalley X to X abs (inelastic, isotropic)
        s = 4;
        v = v_in;
        m = 0.85*m0;

      elseif (rg <= GTot_X(9, E_int))
        % Intervalley X to X em (inelastic, isotropic)
        s = 5;
        v = v_in;
        m = 0.85*m0;

      else
        % Self-scattering -- do nothing
        s = 6;
        v = v_in;
        m = 0.85*m0;
      end

    case 3 % L
      rg = r*GTot_L_Max;
      if (rg <= GTot_L(1, E_int))
        % acoustic (elastic, isotropic)
        s = 1;
        v = v_in;
        m = 0.85*m0;

      elseif (rg <= GTot_L(2, E_int))
        % POP abs (inelastic, anisotropic)
        s = 2;
        v = v_in;
        m = 0.85*m0;

      elseif (rg <= GTot_L(3, E_int))
        % POP em (inelastic, anisotropic)
        s = 3;
        v = v_in;
        m = 0.85*m0;

      elseif (rg <= GTot_L(4, E_int))
        % Intervalley L to G abs (inelastic, isotropic)
        s = 4;
        v = 1;
        m = 0.067*m0;

      elseif (rg <= GTot_L(5, E_int))
        % Intervalley L to G em (inelastic, isotropic)
        s = 5;
        v = 1;
        m = 0.067*m0;

      elseif (rg <= GTot_L(6, E_int))
        % Intervalley L to L abs (inelastic, isotropic)
        s = 4;
        v = v_in;
        m = 0.85*m0;

      elseif (rg <= GTot_L(7, E_int))
        % Intervalley L to L em (inelastic, isotropic)
        s = 5;
        v = v_in;
        m = 0.85*m0;

      elseif (rg <= GTot_L(8, E_int))
        % Intervalley L to X abs (inelastic, isotropic)
        s = 4;
        v = 2;
        m = 0.85*m0;

      elseif (rg <= GTot_L(9, E_int))
        % Intervalley L to X em (inelastic, isotropic)
        s = 5;
        v = 2;
        m = 0.85*m0;

      else
        % Self-scattering -- do nothing
        s = 6;
        v = v_in;
        m = 0.85*m0;
      end

    otherwise
      % do nothing;

  end % case statement


end % function
