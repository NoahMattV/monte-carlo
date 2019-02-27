% ECE 745 Project 2 Part 2
% Noah Van Der Weide

% Rode's Iterative Method
% n-type GaAs
% Include scattering by acoustic phonons, POP, and ionized impurities

% If the scattering mechanism is elastic, then the momentum relaxation rate
% is used. Otherwise, it needs to be iterated. The only mechanism in
% low-field GaAs that isn't elastic is polar-optical-phonon scattering, so
% we use Rode's iterative method to approximate a density of states value
% for the final mobility calculation. 

% The relaxation rates and momentum relaxation rates for the required
% scattering mechanisms were taken from HW 8. 

% For calculation purposes, values were converted to kms values 

% The mobility graph is most likely incorrect, yet I don't know what the
% issue is. 

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

wo = (0.03536*e)/hbar; % longitudinal optical phonon energy hbar*wo (J) (hbar already defined, so explicitly defining wo)
hwo = 0.03536*e; % longitudinal optical phonon energy (J)

No = 1/(exp(hwo/kT) - 1); 
No1 = No + 1;
err = 1e-7;

Dac_G = 7.01*e; % Electron acoustic deformation potential (eV) - Gamma
Dac_L = 9.2*e;
Dac_X = 9.0*e;

mo = 9.11e-31; %kg
m = 0.067*mo;
Ec = 0;
Eg = 1.42*e; % Egap in Joules

Z = 1;
Es = 12.9; % relative dielectric constant
Eo = 8.854e-12; %F/m vacuum permittivity
Ep = Es*Eo; 
%E = linspace(0,2*e,1001); % 0 to 2*e Joules

% To avoid divide-by-zero issues, the energy is set to a very low value in
% lieu of 0. 
E = linspace(1e-255,2*e,1001); % 0 to 2*e Joules
Nd = logspace(14,19,1001); % 10^14 - 10^19 cm^-3
Nd = Nd.*1e6; %convert to m^-3

g = zeros(1001,1001); 
num = zeros(1001);
denom = zeros(1001);

Ipop = zeros(1001);

a = zeros(1001);
b1 = zeros(1001);
b2 = zeros(1001);
fo = exp(-E/kT);

mu = zeros(1001);
t1 = zeros(1001);
t2 = zeros(1001);

% Add up the momentum relaxation rates from the three scattering
% mechanisms for each doping density. 

% ---------------------------------------------------------------------------------------
% 1) Acoustic Phonon Scattering - within elastic and equipartition approximations
% ---------------------------------------------------------------------------------------

gc3d = zeros(1001);
GG = zeros(1001);
GL = zeros(1001);
GX = zeros(1001);
G = zeros(1001);

for i = 1:1001
   %gc3d(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * sqrt(E(i) - Ec); % 3D DOS for GaAs
   %G(i) = 2 * 2*pi/hbar * Dac_G^2 * kT/(2*rho*vs^2) * (1/2) * gc3d(i);
   GG(i) = (Dac_G^2 * kT)/(2*pi*hbar*rho*(vs^2)) * (2*m/(hbar^2))^(3/2) * E(i)^(1/2);
   GL(i) = (Dac_L^2 * kT)/(2*pi*hbar*rho*(vs^2)) * (2*m/(hbar^2))^(3/2) * E(i)^(1/2);
   GX(i) = (Dac_X^2 * kT)/(2*pi*hbar*rho*(vs^2)) * (2*m/(hbar^2))^(3/2) * E(i)^(1/2);
   G(i) = GG(i) + GL(i) + GX(i);
   %G(i) = GG(i);
end

% ---------------------------------------------------------------------------------------
% 2) Polar Optical Phonon Scattering - Absorption and Emission separately
% for electrons in the Gamma valley
% ---------------------------------------------------------------------------------------
% Strongly inelastic (Ge != 0)
% Strongly anisotropic (Gm != G)

Gm_pop_abs = zeros(1001); % momentum relaxation rate
Gm_pop_em = zeros(1001);
Gm_pop = zeros(1001);

Gpop_abs = zeros(1001); % scattering rate
Gpop_em = zeros(1001);
Gpop_tot = zeros(1001);

%wo = (0.03536*e)/hbar; % longitudinal optical phonon energy hbar*wo (eV) (hbar already defined, so explicitly defining wo)
%hwo = 0.03536*e; % longitudinal optical phonon energy (eV)
%No = 1/(exp(hwo/kT) - 1); 

for i = 1:1001   
   Gpop_abs(i) = real(((q^2)*wo*(ko/kinf - 1))/(2*pi*ko*Eo*hbar*sqrt(2*E(i)/m)) * (No*asinh(E(i)/hwo)^(1/2)));
   
   Gm_pop(i) = real(((q^2)*wo*(ko/kinf - 1))/(4*pi*ko*Eo*hbar*sqrt(2*E(i)/m)) * (No*sqrt(1 + hwo/E(i)) + (No+1)*sqrt(1 - hwo/E(i)) - hwo*No/E(i) * asinh(E(i)/hwo)^(1/2) + hwo*(No+1)/E(i) * asinh(E(i)/hwo - 1)^(1/2)));
   
   if (E(i) > hbar*wo)
       Gpop_em(i) = ((q^2)*wo*(ko/kinf - 1))/(2*pi*ko*Eo*hbar*sqrt(2*E(i)/m)) * ((No + 1)*asinh((E(i)/hwo) -1)^(1/2));
   else
       Gpop_em(i) = 0;
   end
   Gpop_tot(i) = Gpop_abs(i) + Gpop_em(i);
end

% ---------------------------------------------------------------------------------------
% 3) Ionized impurity scattering at doping densities of 10^17 and 10^19 cm^-3
% ---------------------------------------------------------------------------------------
% Elastic, Anisotropic

Ld = zeros(1001);
k = zeros(1001);
for i = 1:1001
    k(i) = sqrt(2*m*E(i)/hbar^2);
end

for j = 1:1001
    Ld(j) = sqrt((Ep*kT)/(e^2 * Nd(j))); % Debye Length
end

Gi = zeros(1001,1001);
Gmi = zeros(1001,1001); % momentum relaxation rate

for j = 1:1001 % doping density
    for i = 1:1001 % energy 0-2 eV
        %Gi(i,j) = (Nd(j) * Z^2 * e^4 * Ld(j)^4 * mdos)/(hbar^3 * pi * Eo^2 * Es^2) * ((k(i))/(4*(k(i))^2*Ld(j)^2 + 1));
        
        % Momentum Relaxation Rate
        Gmi(i,j) = (Nd(j)*(Z^2)*(e^4)*m)/(8*pi*(hbar^3)*(Eo^2)*(Es^2)*(k(i)^3)) * (log(1+4*(k(i)^2)*(Ld(j))^2) - (4*(k(i)^2)*(Ld(j)^2))/(1+4*(k(i)^2)*(Ld(j)^2)));    
    end
end

% ---------------------------------------------------------------------------------------
% 4) Adding up all Gammas
% ---------------------------------------------------------------------------------------

Gtot = zeros(1001,1001);

for i = 1:1001 % doping density
    for j = 1:1001 % 0-2 eV
        Gtot(i,j) = G(j) + Gpop_tot(j) + Gmi(i,j);
        %Gtot(i,j) = G(j) + Gm_pop(j) + Gmi(i,j);
    end
end   

% ---------------------------------------------------------------------------------------
% 5) Rode's Iterative Method to calculate g and Ipop
% ---------------------------------------------------------------------------------------
for i = 2:1000 % doping density    
    g(i,1) = 0;
    t1(1) = 0;
    
    % This loop is for defining the density of states values before adding
    % in the Ipop.
    for j = 2:1001
        g(i,j)=(-q/kT*sqrt(2*E(j)/m)*exp(-E(j)/kT))/Gtot(i,j);
        t1(j) = g(i,j);
    end
    
    for j = 2:1000 % Energy
            % The equation for Ipop was found in the supplementary notes
            % breaking the Ipop equation into separate variables
            a(j) = (e^2*wo*m)/(4*pi*Eo*1000*hbar^2*k(j));
            b1(j) = (-1 + (2 - hwo/E(j))/(2*sqrt(1-hwo/E(j))) * log(abs((1+sqrt(1-hwo/E(j)))/(1-sqrt(1-hwo/E(j))))));
            b2(j) = (-1 + (2 + hwo/E(j))/(2*sqrt(1+hwo/E(j))) * log(abs((1+sqrt(1+hwo/E(j)))/(1-sqrt(1+hwo/E(j))))));
        
        % Iteration on g(i,j)
        % Instead of iterating until a minimum error is reached, I opted to
        % iterate ten times. This should give a similar result, as
        % the number of iterations required for g to converge should be
        % small.
        
        for x = 1:10 
            if (E(j) > hwo)
                Ipop(j) = a(i)*(No*g(i,j-1)*b1(j) + No1*g(i,j+1)*b2(j));
            else
                Ipop(j) = a(i)*(No1*g(i,j+1)*b2(j));
            end
            
            t2(j) = g(i,j);
            g(i,j) = t1(j) + Ipop(j)/Gtot(i,j);
        end
    end
    
    % The integral rewritten as a sum becomes very similar to part 1. The
    % delta E is most likely not necessary, but was kept for testing. 
    for j = 1:1000
        num(i) = num(i) + (((E(j)*g(i,j)+E(j+1)*g(i,j+1))/2)*((E(j+1)-E(j))/1001)); % take the average of the current energy*g and the next energy*g
        denom(i) = denom(i) + (((sqrt(E(j))*exp(-E(j)/kT) + sqrt(E(j+1))*exp(-E(j+1)/kT))/2)*((E(j+1)-E(j))/1001));
    end
    mu(i) = -(1/3)*sqrt(2/m)*(num(i)/denom(i));
end

% ---------------------------------------------------------------------------------------
% 6) Plotting the mobility graph
% ---------------------------------------------------------------------------------------

Nd = Nd.*1e-6; % convert back to cm^-3
mu = mu.*1e4; % convert to cm^2/Vs

loglog(Nd(3:1001),mu(3:1001));
grid on;
ylabel('Mobility, \mu_n [cm^2/(Vs)]');
xlabel('Doping Density, N_D [cm^-^3]');
title('Electron Mobility vs. Doping - GaAs at T=300K');

