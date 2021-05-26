% Plot charge and electric field distribution 
clear;
clc;

%% Load data for charge density 
load('bg_charge_n.mat'); % bg_charge_n
bg_charge_n = charge_n;
load('bg_charge_p.mat'); % bg_charge_p
bg_charge_p = charge_p;

load('charge_n.mat'); % charge_n
load('charge_p.mat'); % charge_p

% Look at 2 ps time stamp (1 ps after charge creation)
n = charge_n(:,:,3);
p = charge_p(:,:,3);

%% Subtract background charge density 
delta_n = n - bg_charge_n;
delta_p = p - bg_charge_p;

%% Plot histograms for electrons 
figure();
histogram(delta_n, 500);
axis([0, 7.6e17, 0, 80]);
xlabel('Electron density (cm^{-3})');
ylabel('Number of grid elements');
set(gca, 'FontSize', 16);
grid on;

%% See the min value, max value, and average value of n distribution 
n_min = min(delta_n)
n_max = max(delta_n)
n_mean = mean(delta_n)

%% Plot histograms for holes 
figure();
histogram(delta_p, 500);
axis([0, 1.25e18, 0, 50]);
xlabel('Hole density (cm^{-3})');
ylabel('Number of grid elements');
set(gca, 'FontSize', 16);
grid on;

%% See the min value, max value, and average value of p distribution 
p_min = min(delta_p)
p_max = max(delta_p)
p_mean = mean(delta_p)

%% Input material properties for index calculation - For CdSe 
% All with MKS units 
eps_0 = 8.85e-12; % Vacuum permittivity, in unit of F/m 
eps_m = 2.5^2 * eps_0; % CdTe permittivity CHANGE HERE
n_m = sqrt(eps_m/eps_0); % CdTe refractive index 
e_q = 1.60217662e-19; % Electron charge, in unit of C 
lambda = 1.064e-6; % Target wavelength, in unit of m
c = 3e8; % Speed of light in vacuum, in unit of m/s 
w = 2*pi*c/lambda; % Angular velocity 
mee = 0.13; % Relative electron effective mass Changed 
mhh = 0.45; % Relative hole effective mass  0.4 for CdTe???
m0 = 9.109e-31; % Electron rest mass, in unit of kg 
me = mee*m0; % Electron effective mass 
mh = mhh*m0; % Hole effective mass 
ue = 720e-4; % Electron mobility, with unit of m^2/V*s CHECK
uh = 75e-4; % Hole mobility, with unit of m^2/V*s CHECK

%% Free carrier absorption
free_carrier_abs = 1e-16;
free_carrier_abs_coeff = free_carrier_abs * (mean(delta_n) + mean(delta_p))

%% Free carrier plasma induced refraction
free_carrier_ind = 5.5e-21; % in cm^3
delta_n_fc = free_carrier_ind *  (mean(delta_n) + mean(delta_p))

%% Generalized Drude (Plasma) model 
% This uses the overall carrier density (the perturbed situation) 
% From simulation
Drude_n = mean(n) * 1e6; % Convert to m^-3
Drude_p = mean(p) * 1e6; % Convert to m^-3 

% From analytical
% Drude_n = 1.27e11 * 1e6; % Convert to m^-3
% Drude_p = 1.27e11 * 1e6; % Convert to m^-3 

n_complex = sqrt((eps_m - e_q^2/w*(Drude_n/(me*w + 1i*e_q/ue) + Drude_p/(mh*w + 1i*e_q/uh)))/eps_0);
refraction = real(n_complex);
k_coeff = imag(n_complex);
abs_coeff = 4*pi/(lambda*1e2)*k_coeff
refraction_diff = refraction - sqrt(eps_m/eps_0)

%% Drude Expansion Model
% This uses delta charge density  
% From simulation
Drude_n = mean(delta_n) * 1e6; % Convert to m^-3
Drude_p = mean(delta_p) * 1e6; % Convert to m^-3 

% From analytical
% Drude_n = 1.27e11 * 1e6; % Convert to m^-3
% Drude_p = 1.27e11 * 1e6; % Convert to m^-3 

refraction_diff = -(e_q^2*lambda^2/(8*pi^2*c^2*eps_0*n_m)) * (Drude_n/me + Drude_p/mh)
abs_coeff = -(e_q^3*lambda^2/(4*pi^2*c^3*eps_0*n_m)) * (Drude_n/me^2/ue + Drude_p/mh^2/uh) * 1e-2 % Convert to cm^-1

%% Load data for electric field 
load('bg_electro_E.mat'); % bg_electro_E
bg_electro_E = electro_E;
load('electro_E.mat'); % electro_E

%% Subtract background electric field 
% Look at 2 ps time stamp
Ex = electro_E(:,1,3);
Ey = electro_E(:,2,3);
Ez = electro_E(:,3,3);

bgEx = bg_electro_E(:,1);
bgEy = bg_electro_E(:,2);
bgEz = bg_electro_E(:,3);

delta_Ex = Ex - bgEx;
delta_Ey = Ey - bgEy;
delta_Ez = Ez - bgEz;
%% Plot histograms 
figure();
histogram(delta_Ex, 500);
axis([-3000, 3000, 0, 900]);
xlabel('Electric field strength (V/m)');
ylabel('Number of grid elements');
title('Ex');
set(gca, 'FontSize', 16);
grid on;

%%
figure();
histogram(delta_Ey, 500);
axis([-2500, 2500, 0, 1200]);
xlabel('Electric field strength (V/m)');
ylabel('Number of grid elements');
title('Ey');
set(gca, 'FontSize', 16);
grid on;

%%
figure();
histogram(delta_Ez, 500);
axis([-3.5e4, 3.5e4, 0, 700]);
xlabel('Electric field strength (V/m)');
ylabel('Number of grid elements');
title('Ez');
set(gca, 'FontSize', 16);
grid on;

%% See the min value, max value, and average value of E field distribution 
Ex_min = min(abs(delta_Ex))
Ex_max = max(abs(delta_Ex))
Ex_mean = mean(abs(delta_Ex))

Ey_min = min(abs(delta_Ey))
Ey_max = max(abs(delta_Ey))
Ey_mean = mean(abs(delta_Ey))

Ez_min = min(abs(delta_Ez))
Ez_max = max(abs(delta_Ez))
Ez_mean = mean(abs(delta_Ez))

%% With sign 
Ex_mean = mean(delta_Ex)
Ey_mean = mean(delta_Ey)
Ez_mean = mean(delta_Ez)

%% Pockels effect for CdTe
pockels_coeff = 6.8; % in pm/V
delta_n_pockels = Ez_mean * pockels_coeff * 1e-12

%% FKE
FKE_coeff_alpha = 0.1 / (58e3/1e-2); % E in V/m, for 0.78um
FKE_delta_alpha = Ez_mean * FKE_coeff_alpha % in cm^-1


