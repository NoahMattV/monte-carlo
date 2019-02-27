% ECE 745 Project 2 - Part 1
% Noah Van Der Weide
 
% Relaxation-time approximation (RTA) for low-field transport in Si

% Finding mobility as a function of doping density for Si

% Sweep doping density range from 10^14 - 10^19 cm^-3 and energies 0 to 2eV
% Calculate all Gamma(E) for relevant scattering mechanisms and add them
% together to get Gamma_total(E). Invert Gamma and compute <<T>>. Mobility
% can be found from <<T>>. 

% Gammas for momentum relaxation rate Gm. If the scattering mechanism is
% isotropic (acoustic and intervalley), then Gamma_m = Gamma. 

% For Si, relevant scattering mechanisms in low-field transport are 
% 1) Acoustic Phonon Scattering
% 2) Intervalley Scattering (g and f)
% 3) Ionized Impurity Scattering (which depends on doping concentration!)

% All energies are in Joules
% Units are K/M/S

% The magitude of the final graph is off, but the curve is correct as far
% as I know. 

close all;
clear;

e = 1.602e-19; % charge of an electron/proton
q = e;
hbar = 1.054e-34; % (J*s)
rho = 2329; % (kg/m^3)
kT = 0.0259*e; %(J)
vs = 9040; % longitudinal acoustic velocity (m/s)
Dac = 9.5 * e; % electron acoustic deformation potential (J)
mo = 9.11e-31; % kg
ml = 0.98*mo;
mt = 0.19*mo;
mdos = (0.98*0.19^2)^(1/3)*mo;
Ec = 0;
a = 5.43e-10; %m
Z = 1; 
Es = 11.7; % relative dielectric constant
Eo = 8.854e-12; %F/m vacuum permittivity
Ep = Es*Eo; 

% E = linspace(0,2*e,1001); % 0 to 2*e Joules
E = linspace(1e-255,2*e,1001); % 0 to 2*e Joules

Nd = logspace(14,19,1001); % cm^-3
Nd = Nd.*1e6; %convert to m^-3
Gtot = zeros(54,1001); %doping concentration,energy
Ttot = zeros(54,1001);

mu = zeros(1001);

mc = 0.26*mo;
% ---------------------------------------------------------------------------------------
% Acoustic Phonon Scattering - within elastic and equipartition approximations
% ---------------------------------------------------------------------------------------

gc3d = zeros(1001);
Gac = zeros(1001);

for i = 1:1001
   %gc3d(i) = 6/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * sqrt(E(i) - Ec);
   %G(i) = 2 * 2*pi/hbar * Dac^2 * kT/(2*rho*vs^2) * (1/2) * gc3d(i);
   Gac(i) = (Dac^2 * kT)/(2*pi*hbar*rho*(vs^2)) * (2*mdos/(hbar^2))^(3/2) * E(i)^(1/2);
end

% ---------------------------------------------------------------------------------------
% Equivalent Intervalley Scattering - g and f processes, absorption and emission for each
% ---------------------------------------------------------------------------------------

% f processes - 1 final valley in Si
% g processes - 4 final valleys in Si
Giv_abs_g = zeros(1001);
Giv_em_g = zeros(1001);
Giv_both_g = zeros(1001);

Giv_abs_f = zeros(1001);
Giv_em_f = zeros(1001);
Giv_both_f = zeros(1001);

gc3d_abs_g = zeros(1001);
gc3d_em_g = zeros(1001);
gc3d_abs_f = zeros(1001);
gc3d_em_f = zeros(1001);

% g-type X-X
%Dif_TA_g = 0.5e8; %(eV/cm)
%Dif_TA_g = 5e9; % (eV/m)
Dif_TA_g = 5e9 * e; % (J/m)
%Dif_LA_g = 0.8e8;
%Dif_LA_g = 8e9; % (eV/m)
Dif_LA_g = 8e9 * e; % (J/m)
%Dif_LO_g = 11e8;
%Dif_LO_g = 110e9; %(eV/m)
Dif_LO_g = 110e9 * e; %(J/m)
Dif_g = (Dif_TA_g + Dif_LA_g + Dif_LO_g) / 3;

% f-type X-X
%Dif_TA_f = 0.3e8;
%Dif_TA_f = 3e9; %(eV/m)
Dif_TA_f = 3e9 * e; %(J/m)
%Dif_LA_f = 2e8;
%Dif_LA_f = 20e9; %(eV/m)
Dif_LA_f = 20e9 * e; %(J/m)
%Dif_LO_f = 2e8;
%Dif_LO_f = 20e9; %(eV/m)
Dif_LO_f = 20e9 * e; %(J/m)
Dif_f = (Dif_TA_f + Dif_LA_f + Dif_LO_f) / 3;

Zf_f = 1; % Number of final valleys
Zf_g = 4;

delta_Eif = 0; % 0 for equivalent valley scattering in Si
Eiv_g = e*(0.012 + 0.019 + 0.062) / 3; % (J) %0.062? 
Eiv_f = e*(0.019 + 0.047 + 0.059) / 3;

wif_f = Eiv_f/hbar;
wif_g = Eiv_g/hbar;

Nif_f = 1/(exp(Eiv_f/kT) - 1); % given by Bose-Einstein Distribution Function
Nif_g = 1/(exp(Eiv_g/kT) - 1);

for i = 1:1001
   
   gc3d_abs_f(i) = 6/(2*pi^2) * (2*ml/hbar^2)^(3/2) * real(sqrt(E(i) + Eiv_f - delta_Eif));
   %gc3d_em_f(i) = 6/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - Eiv_f - delta_Eif));
   
   gc3d_abs_g(i) = 6/(2*pi^2) * (2*mt/hbar^2)^(3/2) * real(sqrt(E(i) + Eiv_g - delta_Eif));
   %gc3d_em_g(i) = 6/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - Eiv_g - delta_Eif));
   
   Giv_abs_f(i) = (pi*Dif_f^2*Zf_f)/(2*rho*wif_f) * Nif_f * gc3d_abs_f(i); %absorption
   %Giv_em_f(i) = (pi*Dif_f^2*Zf_f)/(2*rho*wif_f) *(Nif_f + 1) * gc3d_em_f(i); %emission
   
   Giv_abs_g(i) = (pi*Dif_g^2*Zf_g)/(2*rho*wif_g) * Nif_g * gc3d_abs_g(i);
   %Giv_em_g(i) = (pi*Dif_g^2*Zf_g)/(2*rho*wif_g) *(Nif_g + 1) * gc3d_em_g(i);
   
   if (E(i) >= Eiv_f)
       gc3d_em_f(i) = 6/(2*pi^2) * (2*ml/hbar^2)^(3/2) * real(sqrt(E(i) - Eiv_f - delta_Eif));
       Giv_em_f(i) = (pi*Dif_f^2*Zf_f)/(2*rho*wif_f) *(Nif_f + 1) * gc3d_em_f(i); %emission
   else
       Giv_em_f(i) = 0;
   end
   
   if (E(i) >= Eiv_g)
       gc3d_em_g(i) = 6/(2*pi^2) * (2*mt/hbar^2)^(3/2) * real(sqrt(E(i) - Eiv_g - delta_Eif));
       Giv_em_g(i) = (pi*Dif_g^2*Zf_g)/(2*rho*wif_g) *(Nif_g + 1) * gc3d_em_g(i);
   else
       Giv_em_g(i) = 0;
   end
   
   Giv_both_f(i) = Giv_abs_f(i) + Giv_em_f(i);
   Giv_both_g(i) = Giv_abs_g(i) + Giv_em_g(i);
   
end

% ---------------------------------------------------------------------------------------
% Ionized impurity scattering at doping densities of 10^17 and 10^19 cm^-3
% ---------------------------------------------------------------------------------------

k = zeros(1001);
for i = 1:1001
    k(i) = sqrt(2*mdos*E(i)/hbar^2);
end

Ld = zeros(1001);

% calculate the Debye Length for each doping concentration (m)
for i = 1:1001
    Ld(i) = sqrt((Ep*kT)/(q^2 * Nd(i)));
end

Gi = zeros(1001,1001); % For different doping concentrations
Gmi = zeros(1001,1001); % Momentum Relaxation Rate

for i = 1:1001 % doping concentrations
    for j = 1:1001 % 0-2 eV
        %Gi(i,j) = (Nd(i) * Z^2 * q^4 * (Ld(i)^4) * mdos)/(hbar^3 * pi * Eo^2 * Es^2) * (k(j)/(4*(k(j)^2)*(Ld(i)^2) + 1));
        Gmi(i,j) = (Nd(i)*(Z^2)*(e^4)*mdos)/(8*pi*(hbar^3)*(Eo^2)*(Es^2)*(k(j)^3)) * (log(1+4*(k(j)^2)*(Ld(i))^2) - (4*(k(j)^2)*(Ld(i)^2))/(1+4*(k(j)^2)*(Ld(i)^2)));
    
    end
end
%Gmi(1,1) = 0;
% ---------------------------------------------------------------------------------------
% Add up all three scattering mechanisms and invert
% ---------------------------------------------------------------------------------------
T = zeros(1001);
num = zeros(1001);
denom = zeros(1001);

for i = 1:1001 % doping density
    for j = 1:1001 % 0-2 eV
        Gtot(i,j) = Gac(j) + Giv_abs_f(j) + Giv_abs_g(j) + Giv_em_f(j) + Giv_em_g(j) + Gmi(i,j);
        Ttot(i,j) = 1/Gtot(i,j);
        
        % <<T>> = <tau*E>/<E>
        num(i) = num(i) + Ttot(i,j)*E(j)^(3/2)*exp(-E(j)/kT); % <tau*E>
        denom(i) = denom(i) + E(j)^(3/2)*exp(-E(j)/kT); % <E>
        if (E(j) == 0)
            denom(i) = 1;
        end
    end
    T(i) = num(i)/denom(i); %<<Tau>>
    %mobility = e*<<T>>/m* (where m* is the effective conductivity mass -
    %mc = 0.26*mo;
    mu(i) = (q*T(i)/mc); 
end   


mu = mu.*10000; % make cm^2/V-s
Nd = Nd./1e6; % turn back into cm^-3

loglog(Nd(1:1001),mu(1:1001));
ylabel('Mobility, \mu_n [cm^2/(V*s)]');
xlabel('Doping Density, N_D [cm^-^3]');
title('Electron Mobility vs. Doping - Si at T=300K');
grid on;
 
 