% update momentum after scattering for x, y, and z.
function [x_as,y_as,z_as] = getP(x, y, z, p, theta, phi)

  p_as = 1; % I don't know what this is, so I'm setting it to 1. (momentum after scattering)

  % based off of matrix in notes
  sin_theta_cos_phi = sin(theta)*cos(phi); %these may need to change
  sin_theta_sin_phi = sin(theta)*sin(phi);
  cos_theta = cos(theta);

  x_as = (z/p)*(x/sqrt(x^2 + y^2))*sin_theta_cos_phi + (-y/sqrt(x^2 + y^2))*sin_theta_sin_phi + (x/p)*cos_theta * p_as;
  y_as = (y/sqrt(x^2 + y^2))*(z/p)*sin_theta_cos_phi + (x/sqrt(x^2 + y^2))*sin_theta_sin_phi + (y/p)*cos_theta * p_as;
  z_as = (-sqrt(x^2 + y^2)/p)*sin_theta_cos_phi + (z/p)*cos_theta * p_as;

end
