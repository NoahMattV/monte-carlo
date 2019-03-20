% update theta

function out = getTheta(scatt_mech, E_old, E_new)

  rtheta = 2*pi*rand; % between 0 and 2pi
  out = arccos(1-2*rtheta);

  f = 0;
  cos_theta = 0;

  if (scatt_mech == POP_ABS || scatt_mech == POP_EM)
    f = 2*sqrt(E_old*E_new)/((sqrt(E_old) - sqrt(E_new))^2)
    cos_theta = (1 + f - (1 - 2*f)^rtheta)/f;
  end

end
