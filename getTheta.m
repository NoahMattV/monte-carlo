% update theta

function out = getTheta(scatt_mech, E_old, E_new)
  rtheta = 2*pi*rand; % between 0 and 2pi
  f = 0;

  if (scatt_mech == 2 || scatt_mech == 3)
    f = 2*sqrt(E_old*E_new)/((sqrt(E_old) - sqrt(E_new))^2)
    out = arccos((1 + f - (1 - 2*f)^rtheta)/f);
  else
    out = arccos(1-2*rtheta);
  end

end
