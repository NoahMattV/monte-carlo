% ECE 845 Transport in Semiconductor Devices
% Noah Van Der Weide, 2/27/2019
%
% GaAs slab
% Scattering Mechs
% Acoustic, Ionized impurity, POP, and non-equivalent intervalley (Gamma, L, X)


close all;
clear;
numOfParticles = 1001;
numOfTimeSteps = 1001;
e = 1.602e-19;
q = e;
ko = 12.9; % low freq. dielectric const.
kinf = 10.92; % high freq. dielectric const.
hbar = 1.054e-34; % Reduced planck's constant (J*s)
rho = 5360; % density (kg/m^3)
a = 5.462e-10; % lattice constant (m)
kT = 0.0259*e; % (J)
vs = 5240; % Longitudinal acoustic velocity (m/s)

% -----------------------
% Initializing parameters
% -----------------------
E = zeros(numOfParticles);
Px = zeros(numOfParticles);
Py = zeros(numOfParticles);
Pz = zeros(numOfParticles);
valley = ones(numOfParticles); % 1 = Gamma, 2 = X, 3 = L (all start in Gamma)
tff = zeros(numOfParticles);

E_avg = zeros(numOfParticles);
Px_avg = zeros(numOfParticles);
Py_avg = zeros(numOfParticles);
Pz_avg = zeros(numOfParticles);
valley_avg = zeros(numOfParticles);

for i = 1:numOfParticles
  theta = getTheta(); % may have to change this
  phi = getPhi(); % may have to change this

  tff(i) = getTff(); % initial free-flight time for each particle
  E(i) = getEnergy(); % initial kinetic energy for each particle
  % initialize P
end

% -----------------------
% Generate time frame
% -----------------------
deltaT = 10e-14; % 10 fs for now. Shoot for between 1-10 fs.
x1 = 0;
x2 = (numofTimeSteps - 1) * deltaT + x1;
%y = linspace(x1,x2,n) generates n points. The spacing between the points is (x2-x1)/(n-1).
timeStep = linspace(x1, x2, numOfTimeSteps);

% -----------------------
% Generate Scattering Charts (3 charts - 1 for each starting valley)
% Gamma will have 8
% L and X will have 10
% -----------------------


% -----------------------
% Loops!
% -----------------------

for i = 1:numOfTimeSteps % time-stepping loop
  E_tot = 0;
  Px_tot = 0;
  Py_tot = 0;
  Pz_tot = 0;
  valley_tot = 0;
  
  for j = 1:numOfParticles
      if (tff(j) > deltaT) % no scattering
        tff(j) = tff(j) - deltaT;
        %update momentum, p(itime + 1) = p(itime) - eE*deltaT
      else % scattering!
        % check for scattering (update Enew)
        % p(new) = p(itime) - eE*tff
        % come back with updated post scattering, P_as, E_as
        tff(j) = getTff(); % get a new tff
      end % if statement

  E_tot = E_tot + E(j);
  Px_tot = Px_tot + Px(j);
  Py_tot = Py_tot + Py(j);
  Pz_tot = Pz_tot + Pz(j);
  valley_tot = valley_tot + valley(j);

  end % j loop

  % get average E, Px, Py, Pz, and valley occupation(add up and divide by numOfParticles)
  E_avg(i) = E_tot/numOfParticles;
  Px_avg(i) = Px_tot/numOfParticles;
  Py_avg(i) = Py_tot/numOfParticles;
  Pz_avg(i) = Pz_tot/numOfParticles;
  valley_avg(i) = valley_tot/numOfParticles;
end % i loop
