% ECE 845 Transport in Semiconductor Devices
% Noah Van Der Weide, 2/27/2019
%
% GaAs slab (undoped)
% Scattering Mechs
% Acoustic, Ionized impurity, POP, and non-equivalent intervalley (Gamma, L, X)

clc;
close all;
clear;
format long;
numOfParticles = 1001; % Number of particles being tested
numOfTimeSteps = 101; % Number of timesteps. Each timestep should be between 1-10 fs (defined in code).

%global Efield;
%Efield = [0.5 1 2 5 8 10]; % kV/cm Efield(1) = 0.5 ... Efield(6) = 10.
Efield = 0.5;
% Depending on the scattering mechanism applied to the electron, a variable can be equated to these
% This is useful in determining the theta for the angle after scattering
global AC; % Acoustic (elastic, isotropic)
global POP_ABS; % Polar-Optical Phonon (inelastic, anisotropic)
global POP_EM;
global IV_ABS; % Intervalley (inelastic, isotropic)
global IV_EM;
global SELF; % Self-Scattering

AC = 1; % Acoustic (elastic, isotropic)
POP_ABS = 2; % Polar-Optical Phonon (inelastic, anisotropic)
POP_EM = 3;
IV_ABS = 4; % Intervalley (inelastic, isotropic)
IV_EM = 5;
SELF = 6; % Self-Scattering

global e; % charge on an electron
global q;
global ko; % low freq. dielectric const.
global kinf; % high freq. dielectric const.
global hbar; % Reduced planck's constant (J*s)
global rho; % density (kg/m^3)
global a; % lattice constant (m)
global kT; % (J)
global vs; % Longitudinal acoustic velocity (m/s)
global hwo; % longitudinal optical phonon energy (J)

e = 1.602e-19; % charge on an electron
q = e;
ko = 12.9; % low freq. dielectric const.
kinf = 10.92; % high freq. dielectric const.
hbar = 1.054e-34; % Reduced planck's constant (J*s)
rho = 5360; % density (kg/m^3)
a = 5.462e-10; % lattice constant (m)
kT = 0.0259*e; % (J)
vs = 5240; % Longitudinal acoustic velocity (m/s)
hwo = 0.03536*e; % longitudinal optical phonon energy (J)

% -----------------------
% Initializing parameters
% -----------------------
E = zeros(numOfTimeSteps, numOfParticles);
Eint = ones(numOfTimeSteps, numOfParticles); % integer value of energy normalized to number of discrete energy steps in scattering mechanism arrays
eff_m = zeros(numOfTimeSteps, numOfParticles); % effective mass
scatt_mech = zeros(numOfTimeSteps, numOfParticles); % latest scattering mech to affect the particle
valley = ones(numOfTimeSteps, numOfParticles); % 1 = Gamma, 2 = X, 3 = L (all particles start in Gamma)
theta = zeros(numOfTimeSteps, numOfParticles);
phi = zeros(numOfTimeSteps, numOfParticles);
Px = zeros(numOfTimeSteps, numOfParticles);
Py = zeros(numOfTimeSteps, numOfParticles);
Pz = zeros(numOfTimeSteps, numOfParticles);
P = zeros(numOfTimeSteps, numOfParticles);

tff = zeros(1,numOfParticles);

E_avg = zeros(numOfTimeSteps,1);
Px_avg = zeros(numOfTimeSteps,1);
Py_avg = zeros(numOfTimeSteps,1);
Pz_avg = zeros(numOfTimeSteps,1);
valley_avg = zeros(numOfTimeSteps,1);

% -----------------------
% Initializations
% -----------------------
genScatt1 = scatteringMechs();
genScatt = generateScatt();

theta_i = 0;
phi_i = 0;

for j = 1:numOfParticles
  eff_m(1,j) = 0.067;
  tff(1,j) = getTff(); % initial free-flight time for each particle
  [E(1,j), Eint(1,j)] = getEnergy(); % initial kinetic energy for each particle
  % initialize P
  theta_i = getTheta(0,0,0);
  phi_i = getPhi();
  P(1,j) = sqrt(2*eff_m(1,j)*E(1,j));
  Px(1,j) = P(1,j)*sin(theta_i)*cos(phi_i);

  theta_i = getTheta(0,0,0);
  phi_i = getPhi();
  Py(1,j) = P(1,j)*sin(theta_i)*sin(phi_i);

  theta_i = getTheta(0,0,0);
  %phi_i = getPhi();
  Pz(1,j) = P(1,j)*cos(theta_i);
end

%{
figure();
hold on;
histogram(real(Px));
title('px hist');
hold off;

figure();
hold on;
histogram(real(Py));
title('py hist');
hold off;

figure();
hold on;
histogram(real(Pz));
title('pz hist');
hold off;

figure();
hold on;
histogram(E);
title('E hist');
hold off;

figure();
hold on;
histogram(tff);
title('tff hist');
hold off;

figure();
hold on;
histogram(valley);
title('Initial Valley Occupation');
hold off;

%}

% -----------------------
% Generate time frame
% -----------------------
deltaT = 10e-14; % 10 fs for now. Shoot for between 1-10 fs.
x1 = 0;
x2 = (numOfTimeSteps - 1) * deltaT + x1;
%y = linspace(x1,x2,n) generates n points. The spacing between the points is (x2-x1)/(n-1).
timeStep = linspace(x1, x2, numOfTimeSteps);

% -----------------------
% Loops!
% -----------------------
for i = 1:(numOfTimeSteps-1) % time-stepping loop 
  clc;
  fprintf('%d/%d\n',i,numOfTimeSteps);
  E_tot = 0;
  Px_tot = 0;
  Py_tot = 0;
  Pz_tot = 0;
  valley_tot = 0;
  E_old = 0;

  for j = 1:(numOfParticles)
      if (tff(1,j) > deltaT) % no scattering before next timestep.
        tff(1,j) = tff(1,j) - deltaT;
        E(i+1,j) = E(i,j);
        Eint(i+1,j) = Eint(i,j);
        %update momentum, p(itime + 1) = p(itime) - eE*deltaT
        P(i+1,j) = P(i,j) - e*Efield*deltaT;
        Px(i+1,j) = P(i+1,j)*Px(i,j);
        Py(i+1,j) = P(i+1,j)*Py(i,j);
        Pz(i+1,j) = P(i+1,j)*Pz(i,j);

      else
        % scattering before next timestep!
        % check for scattering (update Enew)
        %E_old = E(i,j);

        [scatt_mech(i,j), valley(i,j), eff_m(i,j)] = getScattering(valley(i,j), Eint(i,j));
        [E(i+1,j), Eint(i+1,j)] = updateEnergy(scatt_mech(i,j), E(i,j));

        theta(i,j) = getTheta(scatt_mech(i,j), E(i,j), E(i+1,j));
        phi(i,j) = getPhi();

        % Update Momentum
        % p(new) = p(itime) - eE*tff

        [Px(i+1,j), Py(i+1,j), Pz(i+1,j)] = getP(Px(i,j), Py(i,j), Pz(i,j), P(i,j), theta(i,j), phi(i,j));
        P(i+1,j) = P(i,j) - e*Efield*tff(1,j); % momentum after scattering
        Px(i+1,j) = P(i+1,j)*Px(i+1,j);
        Py(i+1,j) = P(i+1,j)*Py(i+1,j);
        Pz(i+1,j) = P(i+1,j)*Pz(i+1,j);
        % come back with updated post scattering, P_as, E_as

        % Update Free-Flight Time
        tff(1,j) = getTff(); % get a new free-flight time for particle j

      end % if statement

      E_tot = E_tot + E(i,j);
      Px_tot = Px_tot + Px(i,j);
      Py_tot = Py_tot + Py(i,j);
      Pz_tot = Pz_tot + Pz(i,j);
      valley_tot = valley_tot + valley(i,j);

  end % j loop

  % get average E, Px, Py, Pz, and valley occupation(add up and divide by numOfParticles)
  E_avg(i,1) = E_tot/numOfParticles;
  Px_avg(i,1) = Px_tot/numOfParticles;
  Py_avg(i,1) = Py_tot/numOfParticles;
  Pz_avg(i,1) = Pz_tot/numOfParticles;
  valley_avg(i,1) = valley_tot/numOfParticles;
end % i loop

figure();
hold on;
histogram(valley);
title('Valley Occupation');
hold off;

figure();
hold on;
plot(timeStep,E_avg);
title('Ek Avg');
hold off;

figure();
hold on;
plot(timeStep,abs(Px_avg));
title('Px Avg');
hold off;

figure();
hold on;
plot(timeStep,abs(Py_avg));
title('Py Avg');
hold off;

figure();
hold on;
plot(timeStep,abs(Pz_avg));
title('Pz Avg');
hold off;
