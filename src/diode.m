clear; clc;

%%
%% Constants & parameters definition
%%

T = 300;          % Temperature in Kelvin
k = 8.617e-5;     % Boltzmann constant (eV/K)
e0 = 8.85e-14;    % permittivity of free space (F/cm)
q = 1.602e-19;    % charge on an electron (coul)
KS = 11.8;        % Dielectric constant of Si
ni = 1e10;        % intrinsic carrier conc in Silicon at 300K (cm^-3)
EG= 1.12;         % Silicon band gap (eV)

NA = 1e17;        % total acceptors concetration
ND = 1e14;        % total donors conetration

%%
%% Vbi, xn, xp, W, Emax, LDi, LDn, LDp computation in equilibrium conditions
%%

Vbi = (k*T)*log((NA*ND)/(ni^2))                % volts

xn  = sqrt((2*KS*e0/q)*(NA/(ND*(NA+ND)))*Vbi)  % cm
xp  = (ND/NA)*xn                               % cm

W = xn + xp                                    % cm

Emax = ((q*ND)/(KS*e0))*xn                     % volts/cm