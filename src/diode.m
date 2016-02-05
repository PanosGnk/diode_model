%%        DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEER
%%                UNIVERSITY OF THESSALY
%%
%%           ECE656: PRINCIPLES OF SOLID STATE PHYSICS
%%
%%            MATLAB PROJECT: pn-JUNCTION SIMULATION
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
ni = 1.5e10;      % Intrinsic carrier concentration of Silicon at 300K (cm^-3)
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

L = 20*W;                        % total area length (20 times more than depletion width)
dx = min(Ldn, Ldp)/20;           % differential distance (20 times less than minimum Debye length)
N = round(L/dx);                 % total number of points for simulation
M = ceil(N/2);

% Dop = zeros(1, N);               % doping accros the device
Dop(1:M) = -NA/ni;             % 
Dop(M+1:N) = ND/ni;  %

n0 = zeros(1, N);
p0 = zeros(1, N);
V0 = zeros(1, N);

n0(1:M) = 0.5 .* Dop(1:M) .* (1 - sqrt(1 + 1./((Dop(1:M).^2)./4)));
n0(M+1:N) = 0.5 .* Dop(M+1:N) .* (1 + sqrt(1 + 1./((Dop(M+1:N).^2)./4)));

p0 = 1./n0;
V0 = log(n0);


a(1:N) = 1/dx^2;
b(1:N) = - (2/dx^2 + exp(V0) + exp(-V0));
c(1:N) = 1 / dx^2;
V(1:N) = exp(V0) - exp(-V0) - Dop -V0.*(exp(V0)+exp(-V0));

a(1) = 0; a(N) = 0;
c(1) = 0; c(N) = 0;
b(1) = 1; b(N) = 1;
V(1) = V0(1); V(N) = V0(N);

disp('A) Solving Poisson equation in equilibrium conditions.');
tic;

A = gallery('tridiag',a(2:N),b,c(1:N-1));

flag_conv = 0; % convergence of the Poisson loop
delta_acc = 1E-5;
k_iter= 0;

while(~flag_conv) 
 k_iter = k_iter + 1; 
 fi_old=V0;
 V0 = A\V';
 V0 = V0';
 delta = V0-fi_old;
 
 delta_max=max(abs(delta));
 if(delta_max < delta_acc)
       flag_conv = 1;
 else
    for i = 2: N-1
         b(i) = -(2/dx^2 + exp(V0(i)) + exp(-V0(i)));
         V(i) = exp(V0(i)) - exp(-V0(i)) - Dop(i) - V0(i)*(exp(V0(i)) + exp(-V0(i)));
    end  
    A=gallery('tridiag',a(2:N),b,c(1:N-1));
    end
end

toc;



























