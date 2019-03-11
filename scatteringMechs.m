% ECE 845 HW 1 Prob 4
% Noah Van Der Weide
% RTA for a sample of GaAs 10^17 cm^3 n-type Doping
% Room temperature (T = 300K)

% part b) calculate the mobility using RTA for only
% i) acoustic phonon scattering
% ii) ionized impurity scattering
% iii) POP scattering

% part c) Use Matthiessen's rule to calulate mobility and then conductivity.
% How does it compare to experimental results?

% part d) Do not use Matthiessen’s rule, but rather sum up the
% momentum relaxation rates for every given energy, then invert the sum to get
% the total momentum relaxation rate, and then find the average of that
% time, <<τm>>. With it, calculate the mobility and conductivity. How does it
% compare to that obtained using Matthiessen’s rule? How does it compare to
% experiment/your calculation using Rode’s method?

close all;
clear;
global e = 1.602e-19;
global q = e;
global ko = 12.9; % low freq. dielectric const.
global kinf = 10.92; % high freq. dielectric const.
global hbar = 1.054e-34; % Reduced planck's constant (J*s)
global rho = 5360; % density (kg/m^3)
global a = 5.462e-10; % lattice constant (m)
global kT = 0.0259*e; % (J)
global vs = 5240; % Longitudinal acoustic velocity (m/s)

global wo = (0.03536*e)/hbar; % longitudinal optical phonon energy hbar*wo (J) (hbar already defined, so explicitly defining wo)
global hwo = 0.03536*e; % longitudinal optical phonon energy (J)

global No = 1/(exp(hwo/kT) - 1);
global No1 = No + 1;
err = 1e-7;

global Dac_G = 7.01*e; % Electron acoustic deformation potential (eV) - Gamma
global Dac_L = 9.2*e;
global Dac_X = 9.0*e;

global mo = 9.11e-31; %kg
global m = 0.067*mo;
global Ec = 0;
global Eg = 1.42*e; % Egap in Joules

global Z = 1;
global Es = 12.9; % relative dielectric constant
global Eo = 8.854e-12; %F/m vacuum permittivity
global Ep = Es*Eo;
%E = linspace(0,2*e,1001); % 0 to 2*e Joules

% To avoid divide-by-zero issues, the energy is set to a very low value in
% lieu of 0.
global E = linspace(1e-255,2*e,1001); % 0 to 2*e Joules
global Nd = 1e17; % cm^-3
global Nd = Nd.*1e6; %convert to m^-3

global g = zeros(1001);
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

% ---------------------------------------------------------------------------------------
% i) Acoustic Phonon Scattering - within elastic and equipartition approximations
% ---------------------------------------------------------------------------------------

global gc3d = zeros(1001);
global GG = zeros(1001);
global GL = zeros(1001);
global GX = zeros(1001);
global G = zeros(1001);

for i = 1:1001
   %gc3d(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * sqrt(E(i) - Ec); % 3D DOS for GaAs
   %G(i) = 2 * 2*pi/hbar * Dac_G^2 * kT/(2*rho*vs^2) * (1/2) * gc3d(i);
   GG(i) = (Dac_G^2 * kT)/(2*pi*hbar*rho*(vs^2)) * (2*m/(hbar^2))^(3/2) * E(i)^(1/2);
   GL(i) = (Dac_L^2 * kT)/(2*pi*hbar*rho*(vs^2)) * (2*m/(hbar^2))^(3/2) * E(i)^(1/2);
   GX(i) = (Dac_X^2 * kT)/(2*pi*hbar*rho*(vs^2)) * (2*m/(hbar^2))^(3/2) * E(i)^(1/2);
   G(i) = GG(i) + GL(i) + GX(i);
   %G(i) = GG(i);
end

% Mobility calculation for acoustic phonon Scattering

Tac = zeros(1001);
muac = zeros(1001);

for i = 1:1001 % 0-2 eV
    Tac(i) = 1/G(i);

    % <<T>> = <tau*E>/<E>
    num(i) = num(i) + Tac(i)*E(i)^(3/2)*exp(-E(i)/kT); % <tau*E>
    denom(i) = denom(i) + E(i)^(3/2)*exp(-E(i)/kT); % <E>
    if (E(i) == 0)
        denom(i) = 1;
    end
end

for i = 1:1001
  Tac(i) = num(i)/denom(i); %<<Tau>>
  %mobility = e*<<T>>/m* (where m* is the effective conductivity mass -
  muac(i) = (q*Tac(i)/m);
end

% ---------------------------------------------------------------------------------------
% ii) Ionized impurity scattering at doping densities of 10^17
% ---------------------------------------------------------------------------------------
% Elastic, Anisotropic


k = zeros(1001);
for i = 1:1001
    k(i) = sqrt(2*m*E(i)/hbar^2);
end

Ld = sqrt((Ep*kT)/(e^2 * Nd)); % Debye Length

Gi = zeros(1001);
global Gmi = zeros(1001); % momentum relaxation rate

    for i = 1:1001 % energy 0-2 eV
        %Gi(i,j) = (Nd(j) * Z^2 * e^4 * Ld(j)^4 * mdos)/(hbar^3 * pi * Eo^2 * Es^2) * ((k(i))/(4*(k(i))^2*Ld(j)^2 + 1));

        % Momentum Relaxation Rate
        Gmi(i) = (Nd*(Z^2)*(e^4)*m)/(8*pi*(hbar^3)*(Eo^2)*(Es^2)*(k(i)^3)) * (log(1+4*(k(i)^2)*(Ld^2) - (4*(k(i)^2)*(Ld^2))/(1+4*(k(i)^2)*(Ld^2))));

    end

% Mobility calculation for ionized Impurity

Ti = zeros(1001);
mui = zeros(1001);
num = zeros(1001);
denom = zeros(1001);
for i = 1:1001 % 0-2 eV
    Ti(i) = 1/Gmi(i);

    % <<T>> = <tau*E>/<E>
    num(i) = num(i) + Ti(i)*E(i)^(3/2)*exp(-E(i)/kT); % <tau*E>
    denom(i) = denom(i) + E(i)^(3/2)*exp(-E(i)/kT); % <E>
    if (E(i) == 0)
        denom(i) = 1;
    end
end

for i = 1:1001
    Ti(i) = num(i)/denom(i); %<<Tau>>
    %mobility = e*<<T>>/m* (where m* is the effective conductivity mass -
    mui(i) = (q*Ti(i)/m);
end

% ---------------------------------------------------------------------------------------
% iii) Polar Optical Phonon Scattering - Absorption and Emission separately
% for electrons in the Gamma valley
% ---------------------------------------------------------------------------------------
% Strongly inelastic (Ge != 0)
% Strongly anisotropic (Gm != G)

global Gm_pop_abs = zeros(1001); % momentum relaxation rate
global Gm_pop_em = zeros(1001);
Gm_pop = zeros(1001);

global Gpop_abs = zeros(1001); % scattering rate
global Gpop_em = zeros(1001);
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

% mobility calculation for POP
Tpop = zeros(1001);
mupop = zeros(1001);
num = zeros(1001);
denom = zeros(1001);
for i = 1:1001 % 0-2 eV
    Tpop(i) = 1/Gm_pop(i);

    % <<T>> = <tau*E>/<E>
    num(i) = num(i) + Tac(i)*E(i)^(3/2)*exp(-E(i)/kT); % <tau*E>
    denom(i) = denom(i) + E(i)^(3/2)*exp(-E(i)/kT); % <E>
    if (E(i) == 0)
        denom(i) = 1;
    end
end

for i = 1:1001
  Tpop(i) = num(i)/denom(i); %<<Tau>>
  %mobility = e*<<T>>/m* (where m* is the effective conductivity mass -
  mupop(i) = (q*Tpop(i)/m);
end

% ---------------------------------------------------------------------------------------
% Intervalley Scattering - absorption and emission for each
% ---------------------------------------------------------------------------------------
% G, 1 valley
% L, 4 valleys (1/2 * 4 sides)
% X, 3 valleys (1/2 * 6 sides)

%Dif_XL = 2e8; % (eV/cm)
Dif_XL = 2e10*e;% (J/m)
Ein_XL = 0.058*e; % 0.055, 0.041, 0.017 (J)

% Intervalley Deformation Potential
% D_GL = 10e8; % (eV/cm)
% D_GX = 10e8;
% D_LL = 10e8;
% D_LX = 5e8;
% D_XX = 7e8;
D_GL = 10e10*e; % (J/m)
D_GX = 10e10*e;
D_LL = 10e10*e;
D_LX = 5e10*e;
D_XX = 7e10*e;

% Intervalley Phonon Energy
E_GL = 0.0278*e; % (J)
E_GX = 0.0299*e;
E_LL = 0.0290*e;
E_LX = 0.0293*e;
E_XX = 0.0299*e;

% omegas
w_GL = E_GL/hbar;
w_GX = E_GX/hbar;
w_LL = E_LL/hbar;
w_LX = E_LX/hbar;
w_XX = E_XX/hbar;

% Energy separation between valleys
delta_E_LG = 0.29*e; % (J)
delta_E_XG = 0.45*e; % Gamma to X, but I'm afraid to change it now
delta_E_XL = delta_E_LG - delta_E_XG;
% Number of final valleys
Zf_GL = 4;
Zf_GX = 3;
Zf_LL = 3; % because you can only scatter off of 3
Zf_LX = 3;
Zf_XX = 2; % because you can only scatter off of 2

Zf_LG = 1;
Zf_XG = 1;
Zf_XL = 4;

% Bose-Einstein Dist.
No_GL = 1/(exp(E_GL/kT) - 1);
No_GX = 1/(exp(E_GX/kT) - 1);
No_LL = 1/(exp(E_LL/kT) - 1);
No_LX = 1/(exp(E_LX/kT) - 1);
No_XX = 1/(exp(E_XX/kT) - 1);

% Primitives
gc3d_abs_GL = zeros(1001);
gc3d_em_GL = zeros(1001);
global Giv_abs_GL = zeros(1001);
global Giv_em_GL = zeros(1001);
global Giv_GL = zeros(1001);

gc3d_abs_LG = zeros(1001);
gc3d_em_LG = zeros(1001);
global Giv_abs_LG = zeros(1001);
global Giv_em_LG = zeros(1001);
global Giv_LG = zeros(1001);

gc3d_abs_GX = zeros(1001);
gc3d_em_GX = zeros(1001);
global Giv_abs_GX = zeros(1001);
global Giv_em_GX = zeros(1001);
global Giv_GX = zeros(1001);

gc3d_abs_XG = zeros(1001);
gc3d_em_XG = zeros(1001);
global Giv_abs_XG = zeros(1001);
global Giv_em_XG = zeros(1001);
global Giv_XG = zeros(1001);

gc3d_abs_LL = zeros(1001);
gc3d_em_LL = zeros(1001);
global Giv_abs_LL = zeros(1001);
global Giv_em_LL = zeros(1001);
global Giv_LL = zeros(1001);

gc3d_abs_LX = zeros(1001);
gc3d_em_LX = zeros(1001);
global Giv_abs_LX = zeros(1001);
global Giv_em_LX = zeros(1001);
global Giv_LX = zeros(1001);

gc3d_abs_XL = zeros(1001);
gc3d_em_XL = zeros(1001);
global Giv_abs_XL = zeros(1001);
global Giv_em_XL = zeros(1001);
global Giv_XL = zeros(1001);

gc3d_abs_XX = zeros(1001);
gc3d_em_XX = zeros(1001);
global Giv_abs_XX = zeros(1001);
global Giv_em_XX = zeros(1001);
global Giv_XX = zeros(1001);

for i = 1:1001
   gc3d_abs_GL(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) + E_GL - delta_E_LG));
   gc3d_em_GL(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - E_GL - delta_E_LG));

   gc3d_abs_LG(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) + E_GL + delta_E_LG));
   gc3d_em_LG(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - E_GL + delta_E_LG));

   gc3d_abs_GX(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) + E_GX - delta_E_XG));
   gc3d_em_GX(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - E_GX - delta_E_XG));

   gc3d_abs_XG(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) + E_GX + delta_E_XG));
   gc3d_em_XG(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - E_GX + delta_E_XG));

   gc3d_abs_LL(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) + E_LL));
   gc3d_em_LL(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - E_LL));

   gc3d_abs_LX(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) + E_LX - delta_E_XL));
   gc3d_em_LX(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - E_LX - delta_E_XL));

   gc3d_abs_XL(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - E_LX + delta_E_XL));
   gc3d_em_XL(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) + E_LX + delta_E_XL));

   gc3d_abs_XX(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) + E_XX));
   gc3d_em_XX(i) = 1/(2*pi^2) * (2*mdos/hbar^2)^(3/2) * real(sqrt(E(i) - E_XX));


   Giv_abs_GL(i) = (pi*D_GL^2*Zf_GL)/(2*rho*w_GL) * No_GL * gc3d_abs_GL(i); %absorption
   Giv_em_GL(i) = (pi*D_GL^2*Zf_GL)/(2*rho*w_GL) *(No_GL + 1) * gc3d_em_GL(i); %emission
   Giv_GL(i) = Giv_abs_GL(i) + Giv_em_GL(i);

   Giv_abs_LG(i) = (pi*D_GL^2*Zf_LG)/(2*rho*w_GL) * No_GL * gc3d_abs_LG(i); %absorption
   Giv_em_LG(i) = (pi*D_GL^2*Zf_LG)/(2*rho*w_GL) *(No_GL + 1) * gc3d_em_LG(i); %emission
   Giv_LG(i) = Giv_abs_LG(i) + Giv_em_LG(i);

   Giv_abs_GX(i) = (pi*D_GX^2*Zf_GX)/(2*rho*w_GX) * No_GX * gc3d_abs_GX(i);
   Giv_em_GX(i) = (pi*D_GX^2*Zf_GX)/(2*rho*w_GX) *(No_GX + 1) * gc3d_em_GX(i);
   Giv_GX(i) = Giv_abs_GX(i) + Giv_em_GX(i);

   Giv_abs_XG(i) = (pi*D_GX^2*Zf_XG)/(2*rho*w_GX) * No_GX * gc3d_abs_XG(i);
   Giv_em_XG(i) = (pi*D_GX^2*Zf_XG)/(2*rho*w_GX) *(No_GX + 1) * gc3d_em_XG(i);
   Giv_XG(i) = Giv_abs_XG(i) + Giv_em_XG(i);

   Giv_abs_LL(i) = (pi*D_LL^2*Zf_LL)/(2*rho*w_LL) * No_LL * gc3d_abs_LL(i);
   Giv_em_LL(i) = (pi*D_LL^2*Zf_LL)/(2*rho*w_LL) *(No_LL + 1) * gc3d_em_LL(i);
   Giv_LL(i) = Giv_abs_LL(i) + Giv_em_LL(i);

   Giv_abs_LX(i) = (pi*D_LX^2*Zf_LX)/(2*rho*w_LX) * No_LX * gc3d_abs_LX(i);
   Giv_em_LX(i) = (pi*D_LX^2*Zf_LX)/(2*rho*w_LX) *(No_LX + 1) * gc3d_em_LX(i);
   Giv_LX(i) = Giv_abs_LX(i) + Giv_em_LX(i);

   Giv_abs_XL(i) = (pi*D_LX^2*Zf_XL)/(2*rho*w_LX) * No_LX * gc3d_abs_XL(i);
   Giv_em_XL(i) = (pi*D_LX^2*Zf_XL)/(2*rho*w_LX) *(No_LX + 1) * gc3d_em_XL(i);
   Giv_XL(i) = Giv_abs_XL(i) + Giv_em_XL(i);

   Giv_abs_XX(i) = (pi*D_XX^2*Zf_XX)/(2*rho*w_XX) * No_XX * gc3d_abs_XX(i);
   Giv_em_XX(i) = (pi*D_XX^2*Zf_XX)/(2*rho*w_XX) *(No_XX + 1) * gc3d_em_XX(i);
   Giv_XX(i) = Giv_abs_XX(i) + Giv_em_XX(i);
end

% ---------------------------------------------------------------------------------------
% d) Don't use Matthiessen's rule, just sum up Gm
% ---------------------------------------------------------------------------------------
Gtot = zeros(1001);
Ttot = zeros(1001);
num = zeros(1001);
denom = zeros(1001);


for i = 1:1001 % 0-2 eV
    Gtot(i) = G(i) + Gpop_tot(i) + Gmi(i);
    Ttot(i) = 1/Gtot(i);

    % <<T>> = <tau*E>/<E>
    num(i) = num(i) + Ttot(i)*E(i)^(3/2)*exp(-E(i)/kT); % <tau*E>
    denom(i) = denom(i) + E(i)^(3/2)*exp(-E(i)/kT); % <E>
    if (E(i) == 0)
        denom(i) = 1;
    end
end

for i = 1:1001
  T(i) = num(i)/denom(i); %<<Tau>>
  %mobility = e*<<T>>/m* (where m* is the effective conductivity mass -
  mu(i) = (q*T(i)/m);
end
