% update theta

function out = getTheta(scatt_mech, E_old, E_new)
  rtheta = 2*pi*rand; % between 0 and 2pi
  
  if (scatt_mech == 2 || scatt_mech == 3) % If anisotropic (POP)
    f = 2*sqrt(E_old*E_new)/((sqrt(E_old) - sqrt(E_new))^2);
    out = acos((1 + f - (1 - 2*f)^rtheta)/f);

  else % If isotropic
    out = acos(1-2*rtheta);
  end

end
