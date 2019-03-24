% update energy after scattering based on mechanism

function [energy, energyInt] = updateEnergy(scatt_mech, E_in)
  global hwo;
  global e;
  E = linspace(1e-255,2*e,1001); % 0 to 2*e Joules
  % what is w in this? I'm pretty sure it has to do with applied E field.
  switch scatt_mech
    case 1 % Acoustic (elastic, isotropic)
      energy = E_in;
    case 2 % POP Abs (inelastic, anisotropic)
      energy = E_in + hwo;
    case 3 % POP Em
      energy = E_in - hwo;
    case 4 % IV Abs (inelastic, isotropic)
      energy = E_in + hwo;
    case 5 % IV Em
      energy = E_in - hwo;
    case 6 % self
      energy = E_in;
    otherwise
      % do nothing
  end %switch case

  A = repmat(energy, [1 length(E)]);
  [minValue,energyInt] = min(abs(A-E));

end % function
