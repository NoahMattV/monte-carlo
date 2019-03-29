% ECE 845 Transport in Semiconductor Devices
% Noah Van Der Weide, 2/27/2019
%
% GaAs slab (undoped)
% Scattering Mechs
% Acoustic, Ionized impurity, POP, and non-equivalent intervalley (Gamma, L, X)
% Efield is applied along z-direction
clc;
close all;
clear;
format long;
numOfParticles = 1001; % Number of particles being tested
numOfTimeSteps = 201; % Number of timesteps. Each timestep should be between 1-10 fs (defined in code).
deltaT = 10e-15; % 10 fs for now. Shoot for between 1-10 fs.
%global Efield;

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
m0 = 9.11e-31; % kg
%Efield = [0.5 1 2 5 8 10]; % kV/cm Efield(1) = 0.5 ... Efield(6) = 10.
%Efield = 1000*e;

v_avg_E = zeros(20,1);

EfkV = linspace(1,20,20);
for r = 1:20
  v_tot_E = 0;
  Efield = EfkV(r)*1000;
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

valley_G = zeros(numOfTimeSteps, 1);
valley_X = zeros(numOfTimeSteps, 1);
valley_L = zeros(numOfTimeSteps, 1);

E_avg = zeros(numOfTimeSteps,1);
E_avg_G = zeros(numOfTimeSteps,1);
E_avg_X = zeros(numOfTimeSteps,1);
E_avg_L = zeros(numOfTimeSteps,1);

P_avg = zeros(numOfTimeSteps,1);
Px_avg = zeros(numOfTimeSteps,1);
Py_avg = zeros(numOfTimeSteps,1);
Pz_avg = zeros(numOfTimeSteps,1);

valley_avg = zeros(numOfTimeSteps,1);
v_avg = zeros(numOfTimeSteps,1);

eff_m_avg = zeros(numOfTimeSteps,1);
vx_avg = zeros(numOfTimeSteps,1); % for the x and y components when E field is 5 kV/cm
vy_avg = zeros(numOfTimeSteps,1);

% -----------------------
% Initializations
% -----------------------
genScatt1 = scatteringMechs();
genScatt = generateScatt();

theta_i = 0;
phi_i = 0;

for j = 1:numOfParticles
  eff_m(1,j) = 0.067*m0;
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

% -----------------------
% 1) Plot a histogram of the initial energy distribution and the initial momentum
%    along the z-axis. Assume all electrons start in the Gamma valley.
% -----------------------
P_init = zeros(numOfParticles,1);
E_init = zeros(numOfParticles,1);

for i = 1:numOfParticles
    P_init(i) = Pz(1,i)/e;
    E_init(i) = E(1,i)/e;
end
%{
figure();
hold on;
histogram(abs(P_init));
title('Initial Pz');
xlabel('Momentum (kgm/s)');
ylabel('Number of Electrons');
hold off;

figure();
hold on;
histogram(E_init);
title('Initial Energy');
xlabel('Energy (eV)');
ylabel('Number of Electrons');
hold off;
%}
% -----------------------
% Generate time frame
% -----------------------
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
  fprintf('%d/%d\n',r,20);
  E_tot = 0;
  E_tot_G = 0;
  E_tot_X = 0;
  E_tot_L = 0;
  Px_tot = 0;
  Py_tot = 0;
  Pz_tot = 0;
  P_tot = 0;
  valley_tot = 0;
  E_old = 0;
  eff_m_tot = 0;
  v_tot = 0; % velocity

  for j = 1:(numOfParticles)
      if (tff(1,j) > deltaT) % no scattering before next timestep.
        tff(1,j) = tff(1,j) - deltaT;
        E(i+1,j) = E(i,j);
        Eint(i+1,j) = Eint(i,j);
        %update momentum, p(itime + 1) = p(itime) - eE*deltaT
        P(i+1,j) = P(i,j) - e*Efield*deltaT;
        %P(i+1,j) = P(i,j) - Efield*deltaT;
        Px(i+1,j) = P(i+1,j)*Px(i,j);
        Py(i+1,j) = P(i+1,j)*Py(i,j); %TODO this???
        %Px(i+1,j) = Px(i,j);
        %Py(i+1,j) = Py(i,j);
        Pz(i+1,j) = P(i+1,j)*Pz(i,j);
        eff_m(i+1,j) = eff_m(i,j);
      else
        % scattering before next timestep!
        % check for scattering (update Enew)
        %E_old = E(i,j);

        [scatt_mech(i,j), valley(i+1,j), eff_m(i+1,j)] = getScattering(valley(i,j), Eint(i,j));

        if (scatt_mech(i,j) ~= SELF) %self-scattering
            [E(i+1,j), Eint(i+1,j)] = updateEnergy(scatt_mech(i,j), E(i,j), valley(i,j), valley(i+1,j));
            theta(i,j) = getTheta(scatt_mech(i,j), E(i,j), E(i+1,j));
            phi(i,j) = getPhi();

            % Update Momentum
        	% p(new) = p(itime) - eE*tff

            [Px(i+1,j), Py(i+1,j), Pz(i+1,j)] = getP(Px(i,j), Py(i,j), Pz(i,j), P(i,j), theta(i,j), phi(i,j));
            P(i+1,j) = P(i,j) - e*Efield*tff(1,j); % momentum after scattering
            %P(i+1,j) = P(i,j) - Efield*tff(1,j); % momentum after scattering
            Px(i+1,j) = P(i+1,j)*Px(i+1,j);
            Py(i+1,j) = P(i+1,j)*Py(i+1,j); % TODO this???
            Pz(i+1,j) = P(i+1,j)*Pz(i+1,j);
        % come back with updated post scattering, P_as, E_as
        end

        % Update Free-Flight Time
        tff(1,j) = getTff(); % get a new free-flight time for particle j

      end % if statement

      switch valley(i,j)
        case 1
          E_tot_G = E_tot_G + E(i,j);
        case 2
          E_tot_X = E_tot_X + E(i,j);
        case 3
          E_tot_L = E_tot_L + E(i,j);
        otherwise
          disp('E_tot Error');
      end

      E_tot = E_tot + E(i,j);
      Px_tot = Px_tot + Px(i,j);
      Py_tot = Py_tot + Py(i,j);
      Pz_tot = Pz_tot + Pz(i,j);
      P_tot = P_tot + P(i,j);
      valley_tot = valley_tot + valley(i,j);
      eff_m_tot = eff_m_tot + eff_m(i,j);
      v_tot = v_tot + Pz(i,j)/eff_m(i,j); %%%%
  end % j loop

  % get average E, Px, Py, Pz, and valley occupation(add up and divide by numOfParticles)
  valley_G(i,1) = sum(valley(i,:) == 1);
  valley_X(i,1) = sum(valley(i,:) == 2);
  valley_L(i,1) = sum(valley(i,:) == 3);

  E_avg(i,1) = (E_tot/numOfParticles);
  E_avg_G(i,1) = (E_tot_G/valley_G(i,1));
  E_avg_X(i,1) = (E_tot_X/valley_X(i,1));
  E_avg_L(i,1) = (E_tot_L/valley_L(i,1));

  Px_avg(i,1) = (Px_tot/numOfParticles);
  Py_avg(i,1) = (Py_tot/numOfParticles);
  Pz_avg(i,1) = (Pz_tot/numOfParticles);
  P_avg(i,1) = (P_tot/numOfParticles);
  valley_avg(i,1) = valley_tot/numOfParticles;

  %v_avg(i,1) = - v_tot/numOfParticles;
  %v_avg(i+1,1) = - mean(E(i+1,:) - E(i,:))/mean(e*Efield*tff(1,:));
  %v_avg(i,1) = - mean(E(i+1,:) - E(i,:))/mean(Efield*tff(1,:));
  %v_avg(i,1) = - mean(E(i+1,:) - E(i,:))/(e*Efield*deltaT);
  %v_avg(i+1,1) = - mean(E(i+1,:) - E(i,:))/mean(P(i,:));
  %v_avg(i,1) = abs(Pz_avg(i,1)/eff_m_avg(i,1));


  eff_m_avg(i,1) = eff_m_tot/numOfParticles;

  v_avg(i,1) = real(P_avg(i,1)/eff_m_avg(i,1));
  v_tot_E = v_tot_E + v_avg(i,1);
  vx_avg(i,1) = real(Px_avg(i,1)/eff_m_avg(i,1));
  vy_avg(i,1) = real(Py_avg(i,1)/eff_m_avg(i,1));

end % i loop

% Last value isn't calculated, so to make it easier on myself, I'm just setting
% them to NaN
valley_G(numOfTimeSteps,1) = NaN;
valley_X(numOfTimeSteps,1) = NaN;
valley_L(numOfTimeSteps,1) = NaN;
E_avg(numOfTimeSteps,1) = NaN;
E_avg_G(numOfTimeSteps,1) = NaN;
E_avg_X(numOfTimeSteps,1) = NaN;
E_avg_L(numOfTimeSteps,1) = NaN;
vx_avg(numOfTimeSteps,1) = NaN;
vy_avg(numOfTimeSteps,1) = NaN;
eff_m_avg(numOfTimeSteps,1) = NaN;
% -----------------------
% 2) Plot the time evolution of:
%    a) The average electron velocity along the field direction
%    b) Average electron kinetic energy (for each valley as well as the ensemble as a whole)
%    c) Population of each valley for the uniform electric field of 0.5, 1, 2, 5, 8, and 10 kV/cm
%    d) For the electric field of 5 kV/cm, plot the evolution of the x and y components of the electron
%          velocity as well.
% -----------------------

% <v> = - deltaEk/(e*Efield*tff)
%{
% a) The average electron velocity along the field direction
figure();
hold on;
plot(timeStep,v_avg(:,1));
title('Average Velocity Over Time');
xlabel('Time (s)');
ylabel('Velocity (cm/s)');
hold off;

% b) Average electron kinetic energy (for each valley as well as the ensemble as a whole)
figure();
hold on;
plot(timeStep,E_avg/e);
title('Average Ek Over Time');
plot(timeStep,E_avg_G/e);
plot(timeStep,E_avg_X/e);
plot(timeStep,E_avg_L/e);
legend('Total Avg Ek', '\Gamma', 'X', 'L');
xlabel('Time (s)');
ylabel('E_k (eV)');
hold off;

% c) Population of each valley for the uniform electric field of 0.5, 1, 2, 5, 8, and 10 kV/cm
figure();
hold on;
plot(timeStep,valley_G(:,1));
title('Valley Occupation Over Time')
plot(timeStep,valley_X(:,1));
plot(timeStep,valley_L(:,1));
legend('\Gamma', 'X', 'L');
xlabel('Time (s)');
ylabel('Number of Particles');
hold off;

%{
figure();
hold on;
histogram(valley);
title('Valley Occupation at end');
hold off;
%}
%}
% d) For the electric field of 5 kV/cm, plot the evolution of the x and y
%    components of the electron velocity as well.
if (Efield == 5)
  figure();
  hold on;
  plot(timeStep,vx_avg(:,1));
  plot(timeStep,vy_avg(:,1));
  title('X and Y Components of Avg Velocity');
  legend('v_x', 'v_y');
  xlabel('Time (s)');
  ylabel('Velocity (cm/s)');
  hold off;
end



%{
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
%}

% ---------------------------
% 3) From the time evolution of the drift velocity, kinetic energy in each valley, and valley
%    population, extract characteristic times for the approach of each of these quantities to the
%    steady state (in other words, fit the part of the curve describing approach to the steady state
%    with an exp(− t /exp τ), with τ as the fitting parameter). Plot each of these characteristic times
%    τ as a function of the electric field for the field values as above (0.5, 1, 2, 5, 8, and 10
%    kV/cm). How do the obtained relaxation times for the drift velocity and kinetic energy
%    compare to the appropriate <<τ>>’s, that is, <<τm>> and <<τE>>, respectively?
% ---------------------------




% ---------------------------
% 4) Plot the steady state results for the drift velocity, average electron energy, and the valley
%    population versus the electric field. Vary the field from 0.1 to 10 kV/cm. Extract the value of
%    the low-field mobility for GaAs, and compare with the value found in textbooks (give
%    reference). In your steady state calculations, make sure that all transients have died out before
%    extracting the steady state velocity (this will usually happen after ten or a few tens of
%    picoseconds). In order to combat noise when extracting the steady-state quantities, average
%    each quantity of interest over a few picoseconds once you are comfortably in the steady state.
% ---------------------------

v_avg_E(r,1) = v_tot_E/numOfParticles;


end % EfkV

figure();
plot(EfkV,v_avg_E);




% End of function
