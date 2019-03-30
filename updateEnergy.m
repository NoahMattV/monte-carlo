% update energy after scattering based on mechanism
% valleys: Gamma = 1, X = 2, L = 3
function [energy, energyInt] = updateEnergy(scatt_mech, E_in, old_valley, new_valley)
  global hwo;
  %global w_GL;
  %global w_GX;
  %global w_LL;
  %global w_LX;
  %global w_XX;
  %global hbar;
  
  global E_GL; % (J)
  global E_GX;
  global E_LL;
  global E_LX;
  global E_XX;
  
  global e;
  Ec = [1.42*e 1.90*e 1.71*e]; % minimum energy in Gamma, X, L
  E = linspace(1e-255,2*e,1001); % 0 to 2*e Joules

  switch old_valley
    case 1 % Gamma
      if (new_valley == 1) % to Gamma
        %hwif = hbar*hwo;
        hwif = E_GL;
      elseif (new_valley == 2) % to X
        %hwif = hbar*w_GX;
        hwif = E_GX;
      elseif (new_valley == 3) % to L
        %hwif = hbar*w_GL;
        hwif = E_GL;
      else
        disp('Messed up new_valley');
      end
    case 2 % X
      if (new_valley == 1) % to Gamma
        %hwif = hbar*w_GX;
        hwif = E_GX;
      elseif (new_valley == 2) % to X
        %hwif = hbar*w_XX;
        hwif = E_XX;
      elseif (new_valley == 3) % to L
        %hwif = hbar*w_LX;
        hwif = E_LX;
      else
        disp('Messed up new_valley');
      end
    case 3 % L
      if (new_valley == 1) % to Gamma
        %hwif = hbar*w_GL;
        hwif = E_GL;
      elseif (new_valley == 2) % to X
        %hwif = hbar*w_LX;
        hwif = E_LX;
      elseif (new_valley == 3) % to L
        %hwif = hbar*w_LL;
        hwif = E_LL;
      else
        disp('Messed up new_valley');
      end
    otherwise
      disp('uh oh. Messed up old_valley');
  end % switch

  switch scatt_mech
    case 1 % Acoustic (elastic, isotropic)
        energy = E_in;
    case 2 % POP Abs (inelastic, anisotropic)
        energy = E_in + hwo;
    case 3 % POP Em
        energy = E_in - hwo;
    case 4 % IV Abs (inelastic, isotropic)
        energy = E_in + hwif + Ec(old_valley) - Ec(new_valley);
    case 5 % IV Em
        energy = E_in - hwif + Ec(old_valley) - Ec(new_valley);
    case 6 % self
        energy = E_in;
    otherwise
        disp('no valley?');
      % do nothing
  end %switch case

  A = repmat(energy, [1 length(E)]);
  [minValue,energyInt] = min(abs(A-E));

end % function
