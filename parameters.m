%%% MATLAB:
% This code Solves the Arellano (2008) model using C++ with OpenMP.
 
clc;
clear;

%% Step 0 : we need to compile the MEX code.

mex CXXFLAGS="$CXXFLAGS -fopenmp" LDFLAGS="$LDFLAGS -fopenmp" main.cpp ArellanoEconomy.cpp Utils.cpp

%% Step 1: Assuming you have compiled the C++ code into a MEX file named "main_mex"

% Set up the parameters for the economy:
tic;
params.nthreads = 10;
params.beta = 0.9;
params.ny = 51;
params.nb = 251;
params.b_min = -0.1;
params.b_max = 0.5;
params.r = 1.03;
params.rho = 0.945;
params.sigma = 0.025;
params.theta = 0;
params.m = 3.0;
params.y_bar = 0.9778559;
params.max_iter = 1000;
params.tol = 1e-8;

% Call the MEX function to solve the model
output = main(params);

% Retrieve the output arrays
ygrid = output.Ygrid;
bgrid = output.Bgrid;
p = output.P;
v = output.V;
v_r = output.V_r;
v_d = output.V_d;
q = output.Q;
b_p = output.B_p;
d_p = output.D_p;
toc;

% Transform everything from vectors into matrices (n_y x n_b):
Y = ygrid;
B = bgrid;
P = (reshape(p, [params.ny, params.ny]))';
V = (reshape(v, [params.nb, params.ny]))';
V_r = (reshape(v_r, [params.nb, params.ny]))';
Q = (reshape(q, [params.nb, params.ny]))';
D_p = (reshape(d_p, [params.nb, params.ny]))';
B_p = (reshape(b_p, [params.nb, params.ny]))';
ny = params.ny;
nb = params.nb;

bgrid = linspace(-0.1, 0.5, 251);
plot(bgrid, V(1,:),bgrid, V((ny+1)/2,:),bgrid, V(ny,:));
legend({'low y ', 'mid y', 'high y'},'Location','southwest');
title('Value functions for different endowments');
xlabel('Debt due');
ylabel('V');
set(gcf,'color','w');


%% 
T_sim = 1000; % number of periods for the simulation.
Burn = 200;   % Number of observations not taken into account.

[c, b, y] = simulations(V, V_r, B_p, D_p, Q, T_sim, P, Y, B, ny, nb, Burn);

%display(c);
%display(b);
%display(y);

t = linspace(1, T_sim-Burn, T_sim-Burn);
plot(t, y);
legend('Income','Location','southwest');
title('income Simulation');
xlabel('Debt due');
ylabel('y');
set(gcf,'color','w');

%% plot results:
bgrid = linspace(-0.1, 0.5, 251);
plot(bgrid, V(1,:),bgrid, V((ny+1)/2,:),bgrid, V(ny,:));
legend({'low y ', 'mid y', 'high y'},'Location','southwest');
title('Value functions for different endowments');
xlabel('Debt due');
ylabel('V');
set(gcf,'color','w');




