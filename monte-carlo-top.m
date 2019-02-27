% ECE 845 Transport in Semiconductor Devices
% Noah Van Der Weide, 2/27/2019
%
% GaAs slab
% Scattering Mechs
% Acoustic, Ionized impurity, POP, and non-equivalent intervalley (Gamma, L, X)


close all;
clear;
e = 1.602e-19;
q = e;
ko = 12.9; % low freq. dielectric const.
kinf = 10.92; % high freq. dielectric const.
hbar = 1.054e-34; % Reduced planck's constant (J*s)
rho = 5360; % density (kg/m^3)
a = 5.462e-10; % lattice constant (m)
kT = 0.0259*e; % (J)
vs = 5240; % Longitudinal acoustic velocity (m/s)




%------------------------
% Random numbers (do this for each ensemble, I believe)
%------------------------
re = 2*e*rand; % random electron energy between 0 and 2 eV, converted to Joules
rtheta = rand;
rphi = rand;

theta = arccos(1-2*rtheta);
phi = 2*pi*rphi;
E = -(3/2)*kT*log(re);
