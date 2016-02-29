%%          DEPARTMENT OF ELECTRICAL AND COMPUTER ENGINEER
%%                   UNIVERSITY OF THESSALY
%%
%%           ECE656: PRINCIPLES OF SOLID STATE PHYSICS
%%
%%            MATLAB PROJECT: pn-JUNCTION SIMULATION
%%
%% Instructors: George Panagopoulos (gepanago@gmail.com)
%%              George Stamoulis (georges@inf.uth.gr)
%%
%% Project by:  Nonas Evangelos (vagnonas@gmail.com)              
%%              Zampakika Kleopatra (kleopatrakimath@windowslive.com)
%%              Giannakou Panagiotis-Taxiarchis (gpan75@gmail.com)
%%

clear; clc;

%% 1) Constants & parameters definition
%%

T = 300;          % Temperature (K)
k = 1.38e-23;     % Boltzmann constant (J/K)
q = 1.602e-19;    % Charge on an electron (C)
e0 = 1.05e-12;    % Vacuum permittivity (F/cm)
ni = 1.5e10;      % Intrinsic carrier concentration of Silicon at 300K (cm^-3)

Na = 1e16;        % Total acceptors concetration
Nd = 1e16;        % Total donors conetration

%% 2) Vbi, xn, xp, W, Emax, LDi, LDn, LDp computation in equilibrium conditions
%%

Vt = k*T/q;
Vbi = (k*T/q)*log((Na*Nd)/ni^2);          % "built-in" junction voltage (V)
xn = sqrt((2*e0*Na*Vbi)/(q*Nd*(Na+Nd)));  % n-side width of the pn junction deplition region (cm)
xp = sqrt(2*e0*Nd*Vbi/(q*Na*(Na+Nd)));    % p-side width of the pn junction deplition region (cm)
W = sqrt(2*e0*(Na+Nd)*Vbi/(q*Na*Nd));     % depletion width (cm)
Emax = (q/e0)*xp*Na;                      % maximum value of electric field (volts/cm)
Ldn = sqrt((e0*k*T)/(Nd*q^2));            % n-side extrinsic Debye length
Ldp = sqrt((e0*k*T)/(Na*q^2));            % p-side extrinsic Debye length
Ldi = sqrt((e0*k*T)/(ni*q^2));            % intrinsic Debye length

%% 3) Determine simulation area parameters
%%

L = 10*W;                  % total area length (10 times more than depletion width)
dx = min(Ldn,Ldp)/10;      % differential distance (10 times less than minimum Debye length)
N = round(L/dx);           % total number of points for simulation
M = N/2;                   % half-diode length (to determine doping across the device)
dx = dx/Ldi;               % normalize differential distance according to intrinsic Debye length

%% 4) Solve Poisson equation in equilibrium conditions
%%

% define the variables
Dop = zeros(N,1);      % doping across the semiconductor 
n0 = zeros(N,1);       % initial electrons concentration
p0 = zeros(N,1);       % initial holes concentration

a = zeros(N,1);        % first lower diagonal
b = zeros(N,1);        % first upper diagonal
c = zeros(N,1);        % main diagonal

V = zeros(N,1);        % rhs potential
V0 = zeros(N,1);       % initial potential

% initialize doping across the semiconductor
Dop(1:M) = -Na/ni;
Dop(M+1:N) = Nd/ni;

% inital concentration 
n0(1:M) = 0.5*Dop(1:M).*(1 - sqrt(1+1./((Dop(1:M).^2)/4)));
n0(M+1:N) = 0.5*Dop(M+1:N).*(1 + sqrt(1+1./((Dop(M+1:N).^2)/4)));

% initialize potential, electrons & holes concentration, based on
% requirement for charge neutrality
V0 = log(n0);
n = n0;
p = 1./n0;

% initilize three diagonals of matrix as vectors
a(1:N) = 1/dx^2; a(1) = 0; a(N) = 0;
c(1:N) = 1/dx^2; c(1) = 0; c(N) = 0;
b(1:N) = -(2/dx^2 + (exp(V0)+exp(-V0))); b(1) = 1; b(N) = 1;

% initialize sparse matrix from the three vectors
A = gallery('tridiag', a(2:N), b, c(1:N-1));

% initialize rhs vector
V(1:N) = exp(V0) - exp(-V0) + Dop - V0.*(exp(V0)+exp(-V0));
V(1) = V0(1); V(N) = V0(N);

% initialize parameters for iterative poisson equation solver
V_old = zeros(N,1);
V_new = V0;
itol  = 1e-6;
iter = 0;

% begin the iterative procedure: poisson equation solver
disp('Solving poisson equations in equilibrium conditions');
tic;

while(norm(V_new-V_old)/norm(V_new) >= itol)
    iter = iter +1;
    
    % update rhs vector 
    V(2:N-1) = exp(V_new(2:N-1))-exp(-V_new(2:N-1))-Dop(2:N-1)-V_new(2:N-1).*(exp(V_new(2:N-1))+exp(-V_new(2:N-1)));

    % update matrix
    b(2:N-1) = -(2/dx^2 + exp(V_new(2:N-1)) + exp(-V_new(2:N-1)));
    A = gallery('tridiag', a(2:N), b, c(1:N-1));
    
    % solve for new data
    V_old = V_new;
    V_new = A\V;
end

toc;
fprintf('Convergence after %d steps\n\n', iter);

% results after normalization
x = (1:N).*(dx*Ldi);           % diode points on x-axis
V = V_new.*Vt;                 % potential at each point
E = -diff(V)./(dx*Ldi);        % electirc field at each point
r = -diff(V,2);                % charge density at each point

% plot corresponding figures

disp('Figure 1: Potential V(x)');
figure(1);
plot(x, V);
title('Potential for pn-junction in equilibrium conditions');
xlabel('distance x (cm)');
ylabel('potential V(x) (eV)');

disp('Figure 2: Electric field E(x)');
figure(2);
plot(x(1:N-1), E);
title('Electric field for pn-junction in equilibrium conditions');
xlabel('distance x (cm)');
ylabel('electric field E(x) (V/cm)');

disp('Figure 3: Charge density r(x)');
figure(3);
plot(x(1:N-2), r);
title('Charge density for pn-junction in equilibrium conditions');
xlabel('distance x (cm)');
ylabel('charge density r(x) (Q/cm)');

rerr = 100*norm(Emax-abs(min(E)),2)/norm(Emax,2);
fprintf('\nrelative error of maximum electric field value: %.1f', rerr);

fprintf('\nPress any key to exit...\n');
pause;
