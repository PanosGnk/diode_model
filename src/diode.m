%%        DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEER
%%                UNIVERSITY OF THESSALY
%%
%%           ECE656: PRINCIPLES OF SOLID STATE PHYSICS
%%
%%              MATLAB PROJECT: DIODE SIMULATION
%%
%% Instructors: George Panagopoulos (gepanago@gmail.com)
%%              George Stamoulis (georges@inf.uth.gr)
%%
%% Project by:  Giannakou Panagiotis-Taxiarchis (gpan75@gmail.com)
%%              Nonas Evangelos (vagnonas@gmail.com)
%%              Zampakika Kleopatra (kleopatrakimath@windowslive.com)
%%

clear; clc;

%% 1) Constants & parameters definition
%%

T = 300;          % Temperature (K)
k = 1.38e-23;     % Boltzmann constant (J/K)
e0 = 8.85e-14;    % Vacuum permittivity (F/cm)
q = 1.602e-19;    % Charge on an electron (C)
Ks = 11.8;        % Dielectric constant of Silicon at 300K
ni = 1.0e10;      % Intrinsic carrier concentration of Silicon at 300K (cm^-3)
EG = 1.12;        % Silicon band gap (eV)

NA = 1e16;        % Total acceptors concetration
ND = 1e16;        % Total donors conetration

%% 2) Vbi, xn, xp, W, Emax, LDi, LDn, LDp computation in equilibrium conditions
%%

Vbi = (k*T/q)*log((NA*ND)/(ni^2));              % "built-in" junction voltage (V)
xn = sqrt((2*Ks*e0/q)*(NA/(ND*(NA+ND)))*Vbi);   % n-side width of the pn junction deplition region (cm)
xp = (ND/NA)*xn;                                % p-side width of the pn junction deplition region (cm)
W = xn + xp;                                    % depletion width (cm)
Emax = ((q*ND)/(Ks*e0))*xn;                     % maximum value of electric field (volts/cm)
Ldi = sqrt((Ks*e0*k*T)/(2*ni*q^2));             % intrinsic Debye length
Ldn = sqrt((Ks*e0*k*T)/(2*ND*q^2));             % n-side extrinsic Debye length
Ldp = sqrt((Ks*e0*k*T)/(2*NA*q^2));             % p-side extrinsic Debye length


%% 3) Determine simulation area parameters
%%

L = 20*W;                     % total area length (20 times more than depletion width)
dx = min(Ldn, Ldp)/20;        % differential distance (20 times less than minimum Debye length)
N = round(L/dx);              % total number of points for simulation























