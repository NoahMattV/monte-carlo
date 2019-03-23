% update energy after scattering based on mechanism

function [E, Eint] = updateEnergy(scatt_mech, E_in)
  global hwo;
  global e;
  % what is w in this? I'm pretty sure it has to do with applied E field.
  switch scatt_mech
    case 1 % Acoustic (elastic, isotropic)
      E = E_in;
    case 2 % POP Abs (inelastic, anisotropic)
      E = E_in + hwo;
    case 3 % POP Em
      E = E_in - hwo;
    case 4 % IV Abs (inelastic, isotropic)
      E = E_in + hwo;
    case 5 % IV Em
      E = E_in - hwo;
    case 6 % self
      E = E_in;
    otherwise
      % do nothing
  end %switch case

  Eint = ceil(1001*(E/(2*e)));

end % function
